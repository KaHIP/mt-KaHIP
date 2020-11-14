#pragma once

#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/hash_table.h"

#include <array>

namespace parallel{

#undef PERFORMANCE_STATISTICS

template <typename Key, typename Value, typename Hash,
        bool TGrowable = false, bool Cache = true, size_t SizeFactor = 2>
class AdaptiveHashMap {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;
        using hash_type = typename Hash::hash_type;
        using hash_function_type = Hash;
        using Position = uint32_t;

private:
        using TSelf = AdaptiveHashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>;
        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();

        using strategy_type = typename hash_function_type::strategy_type;

#ifdef PERFORMANCE_STATISTICS
#pragma message("PERFORMANCE_STATISTICS FOR HASH TABLE IS ON")
        static size_t num_access;
        static size_t num_contain;
        static size_t overall_max_size;
        static size_t num_probes;
        static size_t num_find_pos;
#endif
public:
        using Iterator = HashTableIterator<TSelf>;

        friend Iterator;

        static constexpr size_t size_factor = SizeFactor;

        explicit AdaptiveHashMap(const uint64_t max_size, const uint8_t max_key_bits_num) :
                _empty_element(std::numeric_limits<Key>::max()),
                _max_size(std::max(round_up_to_next_power_2(max_size), 16u)),
                _ht_size(_max_size * SizeFactor),
                _ht(_ht_size, std::make_pair(_empty_element, Value())),
                _poses(),
                _hash(0, _ht_size, max_key_bits_num, strategy_type::increase),
                _last_key(_empty_element),
                _last_position(0),
                _max_key_bits_num(max_key_bits_num)
        {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(_max_size);
        }

        static size_t get_max_size_to_fit(uint64_t mem) {
                // (SizeFactor + 1.1) * max_size * sizeof(Element) + sizeof(Position) * max_size Bytes = 16 * 1024 Bytes,
                // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8 + 4 * max_size Bytes
                // is the size of a hash table. We calculate that max_size ~ 560.

                size_t size = round_up_to_previous_power_2( mem / (sizeof(Element) * (size_factor)) );
                return size;
        }

        inline constexpr size_t ht_memory_size() const {
                return _ht.size() * sizeof(Element);
        }

        AdaptiveHashMap(const TSelf&) = default;
        AdaptiveHashMap(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

#ifdef PERFORMANCE_STATISTICS
        ~AdaptiveHashMap() {
                overall_max_size = std::max<uint32_t>(overall_max_size, size());
        }
#endif

        Iterator begin() const {
                return Iterator(*this, 0);
        }

        Iterator end() const {
                return Iterator(*this, size());
        }

        inline uint64_t size() const {
                return _poses.size();
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline Value& operator[](const Key& key) {
#ifdef PERFORMANCE_STATISTICS
                ++num_access;
#endif
                if (TGrowable && size() == _max_size)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = Value();
                        _poses.push_back(pos);
                }

                return _ht[pos].second;
        }

        inline bool contains(const Key& key, Value& value) {
#ifdef PERFORMANCE_STATISTICS
                ++num_contain;
#endif
                size_t pos = findPosition(key);
                Element& elem = _ht[pos];
                if (elem.first != _empty_element) {
                        value = elem.second;
                        return true;
                } else {
                        return false;
                }
        }

        inline bool contains(const Key& key) {
#ifdef PERFORMANCE_STATISTICS
                ++num_contain;
#endif
                return _ht[findPosition(key)].first != _empty_element;
        }

        inline void insert(const Element& elem) {
                insertImpl(elem.first, elem.second);
        }

        inline void insert(const Key& key, const Value& value) {
                insertImpl(key, value);
        }

        inline void clear() {
#ifdef PERFORMANCE_STATISTICS
                overall_max_size = std::max<uint32_t>(overall_max_size, size());
#endif
                for (auto pos : _poses) {
                        _ht[pos].first = _empty_element;
                }

                _poses.clear();

                _last_key = _empty_element;
                _last_position = 0;

                _hash.reset_with_increase_strategy(0, _ht_size);
        }

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _max_size = max_size;
                _ht.resize(_ht_size, std::make_pair(_empty_element, Value()));
                _poses.reserve(_max_size);

                _hash.reset(0, _ht_size);
        }

        inline void swap(TSelf& hash_map) {
                std::swap(_ht_size, hash_map._ht_size);
                std::swap(_max_size, hash_map._max_size);
                _ht.swap(hash_map._ht);
                _poses.swap(hash_map._poses);
                std::swap(_last_key, hash_map._last_key);
                std::swap(_last_position, hash_map._last_position);
                std::swap(_empty_element, hash_map._empty_element);
                std::swap(_hash, hash_map._hash);
        }

#ifdef PERFORMANCE_STATISTICS
        static void reset_statistics() {
                num_access = 0;
                num_contain = 0;
                overall_max_size = 0;
                num_find_pos = 0;
                num_probes = 0;
        }

        static void print_statistics() {
                std::cout << "Num access\t" << num_access << std::endl;
                std::cout << "Num contain\t" << num_contain << std::endl;
                std::cout << "Max overall max size\t" << overall_max_size << std::endl;
                std::cout << "Num find pos\t" << num_find_pos << std::endl;
                std::cout << "Num probes\t" << num_probes << std::endl;
                std::cout << "Average num prob per find pos\t" << (num_probes + 0.0) / num_find_pos << std::endl;
        }

        void print_full_statistics() const {
                std::cout << "Num access\t" << num_access << std::endl;
                std::cout << "Num contain\t" << num_contain << std::endl;
                std::cout << "Max overall max size\t" << overall_max_size << std::endl;
                std::cout << "Size\t" << size() << std::endl;
                std::cout << "Mem (table onle)\t" << _ht.size() * sizeof(Element) << std::endl;
                std::cout << "Num find pos\t" << num_find_pos << std::endl;
                std::cout << "Num probes\t" << num_probes << std::endl;
                std::cout << "Average num prob per find pos\t" << (num_probes + 0.0) / num_find_pos << std::endl;
        }
#endif
private:
        explicit AdaptiveHashMap(const uint64_t max_size, const uint8_t max_key_bits_num, const strategy_type strategy) :
                _empty_element(std::numeric_limits<Key>::max()),
                _max_size(std::max(round_up_to_next_power_2(max_size), 16u)),
                _ht_size(_max_size * SizeFactor),
                _ht(_ht_size, std::make_pair(_empty_element, Value())),
                _poses(),
                _hash(0, _ht_size, max_key_bits_num, strategy),
                _last_key(_empty_element),
                _last_position(0),
                _max_key_bits_num(max_key_bits_num)
        {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(_max_size);
        }

        void resize() {
                TSelf new_hash_map(2 * _max_size, _max_key_bits_num, _hash.get_strategy());

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos].first, _ht[pos].second);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key, const Value& value) {
                if (TGrowable && size() == _max_size)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = value;
                        _poses.push_back(pos);
                }
        }

        void rebuild_with_decrease_strategy() {
                _hash.reset_with_decrease_strategy(0, _ht_size);

                std::vector<Element> all_data;
                all_data.reserve(size());

                for (auto pos : _poses) {
                        all_data.emplace_back(_ht[pos].first, _ht[pos].second);
                        _ht[pos].first = _empty_element;
                }

                for (const auto& data : all_data) {
                        insert(data);
                }
        }

        inline Position findPosition(const Key& key) {
#ifdef PERFORMANCE_STATISTICS
//                if ((num_probes + 0.0) / (num_find_pos + 1) > 10) {
//                        reset_statistics();
//                        rebuild_with_decrease_strategy();
//                        return findPosition(key);
//                }
                ++num_find_pos;
#endif
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                const Position startPosition = _hash(key) & (_ht_size - 1);
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
//                        std::cout << "this key = " << key << ", other = " << _ht[pos].first << std::endl;
//                        std::cout << "this key = " << std::bitset<32>(key) << ", other = " << std::bitset<32>(_ht[pos].first) << std::endl;
//                        std::cout << "this h(key) = " << _hash(key) << ", other = " << _hash(_ht[pos].first) << std::endl;
//                        std::cout << "this h(key) = " <<  std::bitset<32>(_hash(key)) << ", other = " <<  std::bitset<32>(_hash(_ht[pos].first)) << std::endl;
//                        std::cout << "this h(key) & ht = " << (_hash(key) & (_ht_size - 1)) << ", other = " << (_hash(_ht[pos].first) & (_ht_size - 1)) << std::endl;
//                        std::cout << "this h(key) & ht = " <<  std::bitset<32>(_hash(key) & (_ht_size - 1)) << ", other = " <<  std::bitset<32>(_hash(_ht[pos].first) & (_ht_size - 1)) << std::endl;
                }

                for (Position pos = 0; pos < startPosition; ++pos) {
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
                }

                throw std::runtime_error("Hash table overflowed");
        }

        Key _empty_element;
        uint64_t _max_size;
        uint64_t _ht_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        Key _last_key;
        Position _last_position;
        const uint8_t _max_key_bits_num;
};

#ifdef PERFORMANCE_STATISTICS
template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t AdaptiveHashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_access(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t AdaptiveHashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_contain(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t AdaptiveHashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::overall_max_size(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t AdaptiveHashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_probes(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, bool Cache, size_t SizeFactor>
size_t AdaptiveHashMap<Key, Value, Hash, TGrowable, Cache, SizeFactor>::num_find_pos(0);
#endif

//template <typename key_type, typename value_type>
//using adaptive_hash_map = AdaptiveHashMap<key_type, value_type, AdaptiveMurmurHash<key_type>, true, false>;

template <typename key_type, typename value_type>
using adaptive_hash_map = AdaptiveHashMap<key_type, value_type, AdaptiveTabularHash<key_type>, true, false>;

template <typename Key, typename Value, bool TGrowable = false, size_t SizeFactor = 2>
class AdaptiveCuckooHashMap {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;
        using hash_type = uint64_t;
        using Position = uint32_t;

private:
        using TSelf = AdaptiveCuckooHashMap<Key, Value, TGrowable, SizeFactor>;
        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();
        static constexpr size_t size_factor = SizeFactor;

public:
        AdaptiveCuckooHashMap(const uint64_t max_size)
                :       _rnd(0)
                ,       _max_size(max_size)
                ,       _ht_size(_max_size * SizeFactor / bucket_type::bucket_size)
                ,       _ht(_ht_size)
        {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(_ht.size() / SizeFactor);
        }

        static size_t get_max_size_to_fit(uint64_t mem) {
                // mem = 2 * n / 8 * 64 , where 2 * n / 8 is the number of nuckets
                // n = 8 * mem / (2 * 64)
                size_t size = round_up_to_previous_power_2( bucket_type::bucket_size * mem / (sizeof(bucket_type) * (size_factor)) );
                return size;
        }

        AdaptiveCuckooHashMap(const TSelf&) = default;
        AdaptiveCuckooHashMap(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

        inline uint64_t size() const {
                return _poses.size();
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline Value& operator[](const Key& key) {
                if (TGrowable && size() == _max_size)
                        resize();

                Position bucket = 0;
                uint32_t slot = 0;

                bool success = findPosition(key, bucket, slot);
                if (!success) {
                        if (slot != bucket_type::dummy_slot) {
                                _ht[bucket][slot].first = key;
                                _ht[bucket][slot].second = Value();
                                _poses.push_back(bucket);
                        } else {
                                success = insert(bucket, {key, Value()}, bucket, slot);
                                // check if was inserted otherwise throw exception
                                if (!success) {
                                        throw std::runtime_error("Can not insert in cuckoo hash table");
                                }
                        }
                }

                return _ht[bucket][slot].second;
        }

        inline bool contains(const Key& key, Value& value) {
                Position bucket;
                uint32_t slot;

                bool success = findPosition(key, bucket, slot);
                if (success) {
                        value = _ht[bucket][slot].second;
                }
                return success;
        }

        inline void swap(TSelf& hash_map) {
                std::swap(_rnd, hash_map._rnd);
                std::swap(_ht_size, hash_map._ht_size);
                std::swap(_max_size, hash_map._max_size);
                _ht.swap(hash_map._ht);
                _poses.swap(hash_map._poses);
                std::swap(_hashes, hash_map._hashes);
        }

        void clear() {
                for (auto pos : _poses) {
                        _ht[pos].clear();
                }

                _poses.clear();
        }

        void reserve(const uint64_t max_size) {
                _ht_size = max_size * SizeFactor / bucket_type::bucket_size;
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _max_size = max_size;
                _ht.resize(_ht_size);
                _poses.reserve(_ht.size() / SizeFactor);
        }

        size_t ht_memory_size() const {
                return _ht_size * sizeof(bucket_type);
        }

private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses) {
                        for (const auto& elem : _ht[pos]) {
                                if (elem.first != _empty_element) {
                                        new_hash_map[elem.first] = elem.second;
                                }
                        }
                        _ht[pos].clear();
                }

                swap(new_hash_map);
        }

        inline bool findPosition(const Key& key, Position& bucket, uint32_t& slot) const {
                uint32_t free_slot = bucket_type::dummy_slot;
                uint32_t free_bucket = 0;
                for (uint8_t i = 0; i < _num_hash_func; ++i) {
                        bool key_found = findPosition(key, i, bucket, slot);
                        if (key_found) {
                                return true;
                        } else {
                                // save free slot
                                if (free_slot == bucket_type::dummy_slot) {
                                        free_bucket = bucket;
                                        free_slot = slot;
                                }
                        }
                }

                bucket = free_bucket;
                slot = free_slot;
                return false;
        }

        inline bool findPosition(const Key& key, uint8_t hash_num, Position& bucket, uint32_t& slot) const {
                bucket = get_bucket(key, hash_num);
                return _ht[bucket].find(key, slot);
        }

        bool insert(Position start_bucket, const Element& const_elem, Position& final_bucket, uint32_t& final_slot) {
                const uint32_t min_slot = 0;
                const uint32_t max_slot = bucket_type::bucket_size - 1;

                uint32_t rnd_slot = _rnd.random_number(min_slot, max_slot);

                Position bucket = start_bucket;
                uint32_t max_iter = 100;
                Element elem = const_elem;

                do {
                        // swap current key with key in the bucket
                        if (elem.first == const_elem.first) {
                                final_bucket = bucket;
                                final_slot = rnd_slot;
                        }

                        _ht[bucket][rnd_slot].swap(elem);

                        uint8_t hash_ind = _rnd.random_number<uint8_t>(0, _num_hash_func - 1);
                        bucket = get_bucket(elem.first, hash_ind);

                        uint32_t free_slot = _ht[bucket].free_slot();

                        if (free_slot != bucket_type::dummy_slot) {
                                _ht[bucket][free_slot] = elem;
                                _poses.push_back(bucket);

                                return true;
                        } else {
                                rnd_slot = _rnd.random_number(min_slot, max_slot);
                        }
                }
                while (max_iter-- > 0);

                return false;
        }

        Position get_bucket(const Key& key, uint8_t hash_num) const {
                if (hash_num == 0) {
                        return std::get<0>(_hashes)(key) & (_ht_size - 1);
                }

                if (hash_num == 1) {
                        return std::get<1>(_hashes)(key) & (_ht_size - 1);
                }

                if (hash_num == 2) {
                        return std::get<2>(_hashes)(key) & (_ht_size - 1);
                }
        }

        struct bucket_type {
                static constexpr uint32_t dummy_slot = std::numeric_limits<uint32_t>::max() - 1;
                static constexpr uint32_t bucket_size = g_cache_line_size / sizeof(Element);
                std::array<Element, bucket_size> data;

                using iterator = typename std::array<Element, bucket_size>::iterator;

                bucket_type() {
                        clear();
                }

                void clear() {
                        for (auto& elem : data) {
                                elem = std::make_pair(_empty_element, Value());
                        }
                }

                bool find(const Key& key, uint32_t& slot) const {
                        slot = dummy_slot;
                        for (size_t i = 0; i < data.size(); ++i) {
                                if (data[i].first == key) {
                                        slot = i;
                                        return true;
                                }
                                if (data[i].first == _empty_element && slot == dummy_slot) {
                                        slot = i;
                                }
                        }
                        return false;
                }

                uint32_t free_slot() const {
                        for (size_t i = 0; i < data.size(); ++i) {
                                if (data[i].first == _empty_element) {
                                        return i;
                                }
                        }
                        return dummy_slot;
                }

                inline constexpr size_t size() const {
                        return data.size();
                }

                inline const Element& operator[] (uint32_t i) const {
                        return data[i];
                }

                inline Element& operator[] (uint32_t i) {
                        return data[i];
                }

                inline iterator begin() {
                        return data.begin();
                }

                inline iterator end() {
                        return data.end();
                }
        };

        constexpr static uint8_t _num_hash_func = 3;
        static constexpr Key _empty_element = std::numeric_limits<Key>::max();

        random _rnd;
        size_t _max_size;
        size_t _ht_size;
        std::vector<bucket_type> _ht;
        std::vector<uint32_t> _poses;
        std::tuple<TabularHash<Key, 3, 2, 10, false>,
                   TabularHash<Key, 2, 2, 10, false>,
                   TabularHash<Key, 0, 3, 10, false>> _hashes;

};

template <typename Key, typename Value, bool TGrowable, size_t SizeFactor>
constexpr Key AdaptiveCuckooHashMap<Key, Value, TGrowable, SizeFactor>::_empty_element;

template <typename key_type, typename value_type>
using adaptive_cuckoo_hash_map = AdaptiveCuckooHashMap<key_type, value_type, true>;

static void test_cuckoo_hash() {
        std::cout << "Start test" << std::endl;
        random rnd(0);

        int test = 1000;
        while (test--) {
                size_t to_insert = 131072;
                size_t access_unlikely_success = to_insert / 2;

                std::unordered_map<uint32_t, uint32_t> data;
                adaptive_cuckoo_hash_map <uint32_t, uint32_t> data_1(to_insert);
                std::vector<std::pair<uint32_t, uint32_t>> to_insert_check;

                to_insert_check.reserve(to_insert);

                for (size_t i = 0; i < to_insert; ++i) {
                        uint32_t key = rnd.random_number(0u, std::numeric_limits<uint32_t>::max() - 2);
                        uint32_t val = rnd.random_number();
                        to_insert_check.emplace_back(key, val);
                }

                std::sort(to_insert_check.begin(), to_insert_check.end());
                auto it = std::unique(to_insert_check.begin(), to_insert_check.end(), [](auto& a, auto& b) {
                        return a.first == b.first;
                });
                to_insert_check.resize(it - to_insert_check.begin());
                std::random_shuffle(to_insert_check.begin(), to_insert_check.end());
                for (size_t i = 0; i < to_insert; ++i) {
                        uint32_t key, val;
                        std::tie(key, val) = to_insert_check[i];

                        data[key] = val;
                        data_1[key] = val;
                }

                for (size_t i = 0; i < to_insert; ++i) {
                        uint32_t key, val;
                        std::tie(key, val) = to_insert_check[i];

                        auto val_1 = data[key];
                        auto val_2 = data_1[key];

                        ALWAYS_ASSERT(val == val_1);
                        ALWAYS_ASSERT(val == val_2);
                }

                for (size_t i = 0; i < access_unlikely_success; ++i) {
                        uint32_t key = rnd.random_number();

                        bool res = data.find(key) != data.end();
                        uint32_t val;
                        bool res_1 = data_1.contains(key, val);

                        ALWAYS_ASSERT(res == res_1);
                }

                size_t i = 0, j = 0;
                while (i < access_unlikely_success || j < to_insert - access_unlikely_success) {
                        bool type;
                        uint32_t key;
                        uint32_t val;
                        if (i < access_unlikely_success && j < to_insert - access_unlikely_success) {
                                if (rnd.bit()) {
                                        key = to_insert_check[j].first;
                                        val = to_insert_check[j].second;
                                        ++j;
                                        type = true;
                                } else {
                                        key = rnd.random_number(0u, std::numeric_limits<uint32_t>::max() - 2);
                                        ++i;
                                        type = false;
                                }
                        } else {
                                if (i < access_unlikely_success) {
                                        key = rnd.random_number(0u, std::numeric_limits<uint32_t>::max() - 2);
                                        ++i;
                                        type = false;
                                }

                                if (j < to_insert - access_unlikely_success) {
                                        key = to_insert_check[j].first;
                                        val = to_insert_check[j].second;
                                        ++j;
                                        type = true;
                                }
                        }

                        if (type) {
                                auto val_1 = data[key];
                                auto val_2 = data_1[key];

                                ALWAYS_ASSERT(val == val_1);
                                ALWAYS_ASSERT(val == val_2);
                        } else {
                                bool res = data.find(key) != data.end();
                                uint32_t val;
                                bool res_1 = data_1.contains(key, val);

                                ALWAYS_ASSERT(res == res_1);
                        }
                }
        }
        std::cout << "Test passed!" << std::endl;
}

}