#pragma once

#include <type_traits>

#include "data_structure/parallel/cache.h"
#include "data_structure/parallel/adaptive_hash_table.h"

namespace parallel {


template <typename _key_type, typename _value_type>
class hash_table_map {
public:
        static_assert(std::is_integral<_key_type>::value, "key_type shoud be integral");
        using key_type = _key_type;
        using value_type = _value_type;
        using map_type = hash_map<key_type, value_type>;
        //using map_type = HashMap<key_type, value_type, parallel::xxhash<key_type>, true>;

        explicit hash_table_map(uint64_t, uint64_t mem)
                :       start_size(get_max_size_to_fit(mem))
                ,       map(start_size)
        {}

        inline bool contains(key_type key, value_type& value) {
                return map.contains(key, value);
        }

        inline void clear() {
                map.clear();
                map.reserve(start_size);
        }

        inline value_type& operator[](key_type key) {
                return map[key];
        }

        inline size_t memory_size() const {
                return map.ht_memory_size();
        }

        inline size_t size() const {
                return map.size();
        }

#ifdef PERFORMANCE_STATISTICS
        static void reset_statistics() {
                map_type::reset_statistics();
        }

        static void print_statistics() {
                map_type::print_statistics();
        }

        void print_full_statistics() const {
                map_type::print_full_statistics();
        }
#endif

private:
        //static constexpr uint64_t mem = g_l2_cache_size;
        //static constexpr uint64_t mem = g_l3_cache_size;

        size_t start_size;
        map_type map;

        static constexpr size_t size_factor = map_type::size_factor;

        static size_t get_max_size_to_fit(uint64_t mem) {
                return map_type::get_max_size_to_fit(mem);
        }
};

template <typename _key_type, typename _value_type>
class array_map {
public:
        static_assert(std::is_integral<_key_type>::value, "key_type shoud be integral");
        using key_type = _key_type;
        using value_type = _value_type;

        explicit array_map(uint64_t _max_size, uint64_t)
                :       max_size(_max_size)
        {
                init(max_size, 0);
        }

        inline void init(uint64_t _max_size, uint64_t) {
                map.assign(_max_size, non_initialized);
                accessed.reserve(128);
        }

        inline bool contains(key_type key, value_type& value) {
                value = map[key];
                return value != non_initialized;
        }

        inline void clear() {
                for (auto key : accessed) {
                        map[key] = non_initialized;
                }
                accessed.clear();
        }

        inline value_type& operator[](key_type key) {
                if (map[key] == non_initialized) {
                        accessed.push_back(key);
                }
                return map[key];
        }

        inline typename std::vector<value_type>::iterator begin() {
                return map.begin();
        }

        inline typename std::vector<value_type>::iterator end() {
                return map.end();
        }

        inline size_t memory_size() const {
                return map.size() * sizeof(value_type);
        }

        inline size_t size() const {
                return accessed.size();
        }

private:
        static constexpr value_type non_initialized = std::numeric_limits<value_type>::max();

        std::vector<value_type> map;
        std::vector<key_type> accessed;
        uint64_t max_size;
};

template <typename _key_type, typename _value_type>
constexpr typename array_map<_key_type, _value_type>::value_type array_map<_key_type, _value_type>::non_initialized;

template <
        typename _key_type,
        typename _value_type
>
class cache_aware_map_impl {
public:
        using key_type = _key_type;
        using value_type = _value_type;
        using small_map_type = hash_table_map<key_type, value_type>;
        using large_map_type = array_map<key_type, value_type>;

        static_assert(std::is_same<typename small_map_type::value_type, typename large_map_type::value_type>::value,
                      "value_type shoud be the same");

        static_assert(std::is_integral<key_type>::value, "key_type shoud be integral");

        cache_aware_map_impl(uint64_t _max_size, uint64_t _start_size, uint32_t _l2_cache_size)
                :       small_map(_max_size, _start_size)
                ,       max_size(_max_size, _start_size)
                ,       large_map_mem(max_size * sizeof(typename large_map_type::value_type))
                ,       cur_container(cur_container_type::SMALL)
                ,       l2_cache_size(_l2_cache_size)
        {
                try_swap_containers();
        }

        inline bool contains(key_type key, value_type& value) {
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.contains(key, value);
                } else {
                        return large_map.contains(key, value);
                }
        }

        inline void clear() {
                cur_container = cur_container_type::SMALL;
                small_map.clear();
                large_map.clear();
        }

        inline value_type& operator[](key_type key) {
                if (cur_container == cur_container_type::SMALL) {
                        value_type& ref = small_map[key];

                        // check if growth happened and if yes then check
                        // hash table fits into L2 cache
                        if (!try_swap_containers())
                                return ref;
                }
                return large_map[key];
        }

        inline size_t memory_size() const {
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.memory_size();
                } else {
                        return large_map.size() * sizeof(value_type);
                }
        }

        inline size_t size() const {
                if (cur_container == cur_container_type::SMALL) {
                        return small_map.size();
                } else {
                        return large_map.size();
                }
        }

        inline size_t is_max_size() const {
                return small_map.memory_size() >= l2_cache_size / 2;
        }
private:
        enum class cur_container_type {
                SMALL,
                LARGE
        };

        small_map_type small_map;
        large_map_type large_map;
        uint64_t max_size;
        const uint64_t large_map_mem;
        cur_container_type cur_container;
        uint32_t l2_cache_size;

        inline bool try_swap_containers() {
                if (cur_container == cur_container_type::SMALL) {
                        uint64_t small_map_mem = small_map.memory_size();

                        if (small_map_mem <= l2_cache_size && small_map_mem <= large_map_mem) {
                                return false;
                        }

                        large_map.init(max_size);

                        for (const auto& rec : small_map) {
                                large_map[rec.first] = rec.second;
                        }

                        small_map.clear();
                        cur_container = cur_container_type::LARGE;
                        return true;
                }
                return false;
        }
};

template <typename key_type, typename value_type>
using cache_aware_map = hash_table_map<key_type, value_type>;
}
