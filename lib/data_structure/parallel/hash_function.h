#pragma once

#include <bitset>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#define XXH_PRIVATE_API
#include "data_structure/parallel/lib/xxhash.h"
#include "data_structure/parallel/bits.h"
#include "tools/macros_assertions.h"

namespace parallel {

template <typename T>
struct simple_hash {
        static const size_t significant_digits = 64;
        using hash_type = uint64_t;

        simple_hash()
        {}

        explicit simple_hash(uint32_t)
        {}

        inline hash_type operator()(const T& x) const {
                return x;
        }
};

template <typename T>
struct xxhash {
        static const size_t significant_digits = 64;
        using hash_type = uint64_t;

        explicit xxhash(uint32_t seed = 0)
                : _seed(seed)
        {}

        inline void reset(uint32_t seed) {
                _seed = seed;
        }

        inline hash_type operator()(const T& x) const {
                return XXH64(&x, sizeof(x), _seed);
        }

private:
        uint32_t _seed;
};

template <typename T>
struct xxhash<std::pair<T, T>> {
        using hash_type = uint64_t;

        explicit xxhash(uint32_t seed = 0)
                : _h(seed)
        {}

        inline void reset(uint32_t seed) {
                _h.reset(seed);
        }

        uint64_t operator()(const std::pair<T, T>& x) const {
                return _h(x.first) ^ _h(x.second);
        }

private:
        xxhash<T> _h;
};

template <typename Key>
class MurmurHash {
public:
        static const size_t significant_digits = 64;
        using hash_type = std::uint64_t;

        explicit MurmurHash(uint32_t seed = 0)
                : _seed(seed) {}

        inline void reset(uint32_t seed) {
                _seed = seed;
        }

        inline hash_type operator()(const Key& key) const {
                return hash(reinterpret_cast<const void*>(&key), sizeof(key), _seed);
        }

private:
        uint32_t _seed;

        inline hash_type hash(const void* key, uint32_t len, uint32_t seed) const {
                const uint64_t m = 0xc6a4a7935bd1e995;
                const int r = 47;

                uint64_t h = seed ^(len * m);

                const uint64_t* data = (const uint64_t*) key;
                const uint64_t* end = data + (len / 8);

                while (data != end) {
                        uint64_t k = *data++;

                        k *= m;
                        k ^= k >> r;
                        k *= m;

                        h ^= k;
                        h *= m;
                }

                const unsigned char* data2 = (const unsigned char*) data;

                switch (len & 7) {
                        case 7:
                                h ^= uint64_t(data2[6]) << 48;
                        case 6:
                                h ^= uint64_t(data2[5]) << 40;
                        case 5:
                                h ^= uint64_t(data2[4]) << 32;
                        case 4:
                                h ^= uint64_t(data2[3]) << 24;
                        case 3:
                                h ^= uint64_t(data2[2]) << 16;
                        case 2:
                                h ^= uint64_t(data2[1]) << 8;
                        case 1:
                                h ^= uint64_t(data2[0]);
                                h *= m;
                };

                h ^= h >> r;
                h *= m;
                h ^= h >> r;

                return h;
        }
};

template <typename T>
struct MurmurHash<std::pair<T, T>> {
        using hash_type = uint64_t;

        uint64_t operator()(const std::pair<T, T>& x) const {
                return _h(x.first) ^ _h(x.second);
        }

        inline void reset(uint32_t seed) {
                _h.reset(seed);
        }
private:
        MurmurHash<T> _h;
};

template <typename T, uint8_t _key_bits_num = 9, uint8_t _num_tables = 2, uint8_t _table_bits_num = 8,
        bool add_key = true,
        typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
class TabularHash {
public:
        using hash_type = uint64_t;

        /*
         * int bit representation = [..., num_tables * table_bits_num, key_bits_num]
         *                          ^                                              ^
         *                most significant bit                           least significant bit
         */
private:
        // best = 2
        static constexpr uint8_t num_tables = _num_tables;

        // first key_bits_num least significant bits, best = 9
        static constexpr uint8_t key_bits_num = _key_bits_num;

        // table_bits_num bits AFTER key_bits_num least significant bits, best = 8
        static constexpr uint8_t table_bits_num = _table_bits_num;

        static constexpr uint64_t table_size = 1 << table_bits_num;

        static constexpr uint64_t key_mask = (1 << key_bits_num) - 1;
        static constexpr uint64_t table_mask = (1 << table_bits_num) - 1;

        T table[num_tables * table_size];
public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < num_tables * table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> key_bits_num;
                for (uint8_t i = 0; i < num_tables; ++i) {
                        res ^= table[i * table_size + (shifted & table_mask)];
                        shifted >>= table_bits_num;
                }

                if (add_key) {
                        return res ^ (key & key_mask);
                }
                else {
                        return res;
                }
        }
};

template <typename T>
class TabularHash<T, 3, 2, 10, true> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 3;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res ^ (key & 0x7);
        }
};

template <typename T>
class TabularHash<T, 3, 2, 10, false> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 3;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res;
        }
};

template <typename T>
class TabularHash<T, 4, 2, 10, true> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 4;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res ^ (key & 0xF);
        }
};

template <typename T>
class TabularHash<T, 4, 2, 10, false> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 4;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res;
        }
};

template <typename T>
class TabularHash<T, 5, 2, 10, true> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 5;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res ^ (key & 0x1F);
        }
};

template <typename T>
class TabularHash<T, 5, 2, 10, false> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 5;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res;
        }
};

template <typename T, bool add_key>
class TabularHash<T, 0, 2, 10, add_key> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                res ^= table[key & 0x3FF];
                res ^= table[1024 + ((key >> 10) & 0x3FF)];

                return res;
        }
};

template <typename T>
class TabularHash<T, 2, 2, 10, true> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 2;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res ^ (key & 0x3);
        }
};

template <typename T>
class TabularHash<T, 2, 2, 10, false> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 2048;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                T shifted = key >> 2;
                res ^= table[shifted & 0x3FF];
                res ^= table[1024 + ((shifted >> 10) & 0x3FF)];

                return res;
        }
};

template <typename T, bool add_key>
class TabularHash<T, 0, 3, 10, add_key> {
public:
        using hash_type = uint64_t;

private:
        constexpr static size_t table_size = 3072;

        T table[table_size];

public:
        explicit TabularHash(uint32_t seed = 0) {
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < table_size; ++i) {
                        table[i] = rnd(mt);
                }
        }

        uint64_t operator()(T key) const {
                uint64_t res = 0;
                res ^= table[key & 0x3FF];
                res ^= table[1024 + ((key >> 10) & 0x3FF)];
                res ^= table[2048 + ((key >> 20) & 0x3FF)];

                return res;
        }
};

class BaseAdaptiveHash {
public:
        enum class strategy_type {
                increase,
                decrease
        };

        BaseAdaptiveHash(uint8_t _max_full_key_bits_num, size_t _ht_size, strategy_type _strategy)
                :       max_full_key_bits_num(_max_full_key_bits_num),
                        min_key_bits_num(3),
                        max_key_bits_num(9),
                        key_bits_num(0),
                        strategy(_strategy),
                        ht_size(_ht_size),
                        key_mask(0)
        {}

        strategy_type get_strategy() const {
                return strategy;
        }

protected:
        uint8_t max_full_key_bits_num;
        uint8_t min_key_bits_num;
        uint8_t max_key_bits_num;
        uint8_t key_bits_num;
        strategy_type strategy;

        uint32_t ht_size;
        uint64_t key_mask;
        static constexpr uint64_t ht_min_capacity_tabular_hash = 32768;

        void update_key_bits_num() {
                ALWAYS_ASSERT(ht_min_capacity_tabular_hash != 0);
                ALWAYS_ASSERT(is_power_2(ht_size));
                ALWAYS_ASSERT(ht_size >= ht_min_capacity_tabular_hash);

                if (strategy == strategy_type::increase) {
                        key_bits_num = min_key_bits_num + log2(ht_size / ht_min_capacity_tabular_hash);
                        key_bits_num = std::min(key_bits_num, max_key_bits_num);
                        key_mask = (1 << key_bits_num) - 1;
                }

                if (strategy == strategy_type::decrease) {
                        if (key_bits_num <= min_key_bits_num) {
                                key_bits_num = 0;
                        } else {
                                --key_bits_num;
                        }
                }
        }

};

template <typename Key>
class AdaptiveMurmurHash : public BaseAdaptiveHash {
private:
        using base_type = BaseAdaptiveHash;
public:
        using hash_type = std::uint64_t;

        explicit AdaptiveMurmurHash(uint32_t seed, size_t _ht_size, uint8_t _max_full_key_bits_num,
                                    strategy_type _strategy)
                :       base_type(_max_full_key_bits_num, _ht_size, _strategy),
                        hash(seed)
        {
                reset(seed, _ht_size);
        }

        inline hash_type operator()(const Key& key) const {
                Key shifted_key = key >> key_bits_num;
                return hash(shifted_key) ^ (key & key_mask);
        }

        void reset_with_decrease_strategy(uint32_t seed, size_t _ht_size) {
                strategy = strategy_type::decrease;
                reset(seed, _ht_size);
        }

        void reset_with_increase_strategy(uint32_t seed, size_t _ht_size) {
                strategy = strategy_type::increase;
                reset(seed, _ht_size);
        }

        void reset(uint32_t seed, size_t _ht_size) {
                hash.reset(seed);
                ht_size = _ht_size;

                update_key_bits_num();
                std::cout << "strategy\t" << (strategy == strategy_type::increase ? "increase" : "decrease")
                          << std::endl;
                std::cout << "ht_size\t" << ht_size << std::endl;
                std::cout << "key_bits_num\t" << (int) key_bits_num << std::endl;
        }

private:
        MurmurHash<Key> hash;
        uint64_t key_mask;
};

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
class AdaptiveTabularHash : public BaseAdaptiveHash {
private:
        using base_type = BaseAdaptiveHash;

        static constexpr uint8_t max_table_bits_num = 8;
        /*
         * int bit representation = [..., table_bits_num[1], table_bits_num[0], key_bits_num]
         *                          ^                                              ^
         *                most significant bit                           least significant bit
         */

        uint8_t num_tables;
        std::unique_ptr<T[]> table;
        std::unique_ptr<uint8_t[]> table_bits_num;
        std::unique_ptr<uint64_t[]> table_mask;
        std::unique_ptr<uint32_t[]> table_size_prefix;

public:
        using hash_type = uint64_t;

        AdaptiveTabularHash(uint32_t seed, size_t _ht_size, uint8_t _max_full_key_bits_num, strategy_type _strategy)
                : base_type(_max_full_key_bits_num, _ht_size, _strategy)
        {
                reset(seed, _ht_size);
        }

        void reset_with_decrease_strategy(uint32_t seed, size_t _ht_size) {
                strategy = strategy_type::decrease;
                reset(seed, _ht_size);
        }

        void reset_with_increase_strategy(uint32_t seed, size_t _ht_size) {
                strategy = strategy_type::increase;
                reset(seed, _ht_size);
        }

        void reset(uint32_t seed, size_t _ht_size) {
                ALWAYS_ASSERT(max_full_key_bits_num != 0);

                ht_size = _ht_size;

                update_key_bits_num();

                uint8_t total_bits_num = max_full_key_bits_num;
                ALWAYS_ASSERT(total_bits_num >= key_bits_num);

                uint8_t total_table_bits_num = total_bits_num - key_bits_num;
                num_tables =
                        total_table_bits_num / max_table_bits_num + (total_table_bits_num % max_table_bits_num != 0);
                table_bits_num = std::make_unique<uint8_t[]>(num_tables);
                table_mask = std::make_unique<uint64_t[]>(num_tables);
                table_size_prefix = std::make_unique<uint32_t[]>(num_tables);

                uint32_t total_table_size = 0;
                for (size_t i = 0; i < total_table_bits_num / max_table_bits_num; ++i) {
                        table_bits_num[i] = max_table_bits_num;
                        table_mask[i] = (1 << max_table_bits_num) - 1;
                        table_size_prefix[i] = total_table_size;
                        total_table_size += table_mask[i] + 1;
                }

                uint8_t remainder = total_table_bits_num % max_table_bits_num;
                if (remainder) {
                        table_bits_num[num_tables - 1] = remainder;
                        table_mask[num_tables - 1] = (1 << remainder) - 1;
                        table_size_prefix[num_tables - 1] = total_table_size;
                        total_table_size += table_mask[num_tables - 1] + 1;
                }

                table = std::make_unique<T[]>(total_table_size);
                std::uniform_int_distribution<uint32_t> rnd;
                std::mt19937 mt(seed);

                for (size_t i = 0; i < total_table_size; ++i) {
                        table[i] = rnd(mt);
                }

                std::cout << "Constructed tabular hash function:" << std::endl;
                std::cout << "max_full_key_bits_num\t" << (int)max_full_key_bits_num << std::endl;
                std::cout << "total_table_size\t" << total_table_size << std::endl;
                std::cout << "strategy\t" << (strategy == strategy_type::increase ? "increase" : "decrease")
                          << std::endl;
                std::cout << "Ht size\t" << ht_size << std::endl;
                std::cout << "min_key_bits_num\t" << (int) min_key_bits_num << std::endl;
                std::cout << "max_key_bits_num\t" << (int) max_key_bits_num << std::endl;
                std::cout << "total_bits_num\t" << (int) total_bits_num << std::endl;
                std::cout << "num_tables\t" << (int) num_tables << std::endl;
                std::cout << "key_bits_num\t" << (int) key_bits_num << std::endl;
                std::cout << "key_mask\t" << std::bitset<32>(key_mask) << std::endl;

                size_t shift = key_bits_num;
                for (size_t i = 0; i < num_tables; ++i) {
                        std::cout << "table " << i + 1 << " num bits \t" << (int) table_bits_num[i] << std::endl;
                        std::cout << "table " << i + 1 << " mask \t" << std::bitset<32>(table_mask[i] << shift)
                                  << std::endl;
                        shift += table_bits_num[i];
                }
        }

        inline hash_type operator()(T key) const {
//                uint64_t res = 0;
//                T shifted = key >> key_bits_num;
//                for (uint8_t i = 0; i < num_tables; ++i) {
//                        res ^= table[table_size_prefix[i] + (shifted & table_mask[i])];
//                        shifted >>= table_bits_num[i];
//                }
//                return res ^ (key & key_mask);
                uint64_t res = 0;
                T shifted = key >> key_bits_num;
                res ^= table[shifted & 0xFF];
                res ^= table[256 + ((shifted >> 8) & 0xFF)];
                return res ^ (key & key_mask);
        }
};
}
