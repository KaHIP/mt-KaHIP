#pragma once

#include <algorithm>
#include <bitset>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "data_structure/parallel/bits.h"
#include "data_structure/parallel/cache.h"
#include "data_structure/parallel/hash_function.h"



#include "tools/macros_assertions.h"

namespace parallel {
template <typename hash_table_type>
static constexpr size_t get_max_size_to_fit_l1() {
        // (SizeFactor + 1.1) * max_size * sizeof(Element) Bytes = 16 * 1024 Bytes,
        // where 16 * 1024 Bytes is half of L1 cache and (2 + 1.1) * max_size * 8
        // is the size of a hash table. We calculate that max_size ~ 512.
        return round_up_to_previous_power_2(16 * 1024 / (sizeof(typename hash_table_type::Element) * (hash_table_type::size_factor)));
}

template <typename HashTable>
class HashTableIterator {
private:
        using HashMap = HashTable;
        using Element = typename HashMap::Element;
        using Position = typename HashMap::Position;

public:
        HashTableIterator(const HashMap& hm, const Position offset) :
                _hm(hm),
                _offset(offset)
        {}

        const Element& operator* () {
                return _hm._ht[_hm._poses[_offset]];
        }

        HashTableIterator& operator++ () {
                ++_offset;
                return *this;
        }

        bool operator== (const HashTableIterator& it) {
                return _offset == it._offset;
        }

        bool operator!= (const HashTableIterator& it) {
                return !(*this == it);
        }

private:
        const HashMap& _hm;
        Position _offset;
};

template <typename HashTable>
class HashTableWithEraseIterator {
private:
        using HashMapWithErase = HashTable;
        using Element = typename HashMapWithErase::Element;
        using Position = typename HashMapWithErase::Position;

public:
        HashTableWithEraseIterator(const HashMapWithErase& hm, const Position offset) :
                _hm(hm),
                _offset(offset)
        {
                get_next();
        }

        const Element& operator* () {
                return _hm._ht[_hm._poses[_offset]];
        }

        HashTableWithEraseIterator& operator++ () {
                ++_offset;
                get_next();
                return *this;
        }

        bool operator== (const HashTableWithEraseIterator& it) {
                return _offset == it._offset;
        }

        bool operator!= (const HashTableWithEraseIterator& it) {
                return !(*this == it);
        }

private:
        const HashMapWithErase& _hm;
        Position _offset;

        inline void get_next() {
                while (_offset < _hm._poses.size() && _hm._ht[_hm._poses[_offset]].first == _hm._deleted_element) {
                        ++_offset;
                }
        }
};

#undef PERFORMANCE_STATISTICS

template <typename Key, typename Value, typename Hash = xxhash<Key>,
        bool TGrowable = false, size_t SizeFactor = 2>
class HashMap {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;
        using hash_type = typename Hash::hash_type;
        using hash_function_type = Hash;
        using Position = uint32_t;

private:
        using TSelf = HashMap<Key, Value, Hash, TGrowable, SizeFactor>;
        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();
        static constexpr size_t size_factor = SizeFactor;

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

        explicit HashMap(const uint64_t max_size = 1) :
                _empty_element(std::numeric_limits<Key>::max()),
                _max_size(std::max(round_up_to_next_power_2(max_size), 16u)),
                _ht_size(_max_size * SizeFactor),
                //_ht(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value())),
                _ht(_ht_size + 301, std::make_pair(_empty_element, Value())),
                _poses(),
                _hash()
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

        HashMap(const TSelf&) = default;
        HashMap(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

#ifdef PERFORMANCE_STATISTICS
        ~HashMap() {
                overall_max_size = std::max<uint32_t>(overall_max_size, size());
        }
#endif

        Iterator begin() const {
                return Iterator(*this, 0);
        }

        Iterator end() const {
                return Iterator(*this, size());
        }

        void reserve(const uint64_t max_size) {
                _ht_size = max_size * SizeFactor;
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _max_size = max_size;
                //_ht.resize(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value()));

                // +301 and dummy element are to remove check for the size of _ht in findPosition function
                _ht.resize(_ht_size + 301, std::make_pair(_empty_element, Value()));

                _poses.reserve(_max_size);
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

                if (pos + 1 == _ht.size()) {
                        resize();
                        return (*this)[key];
                }

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
                value = elem.second;
                return elem.first != _empty_element;
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
        }

        inline void swap(TSelf& hash_map) {
                std::swap(_ht_size, hash_map._ht_size);
                std::swap(_max_size, hash_map._max_size);
                _ht.swap(hash_map._ht);
                _poses.swap(hash_map._poses);
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
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos].first, _ht[pos].second);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key, const Value& value) {
                if (TGrowable && size() == _max_size)
                        resize();

                const Position pos = findPosition(key);

                if (pos + 1 == _ht.size()) {
                        resize();
                        insertImpl(key, value);
                        return;
                }

                if (_ht[pos].first == _empty_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = value;
                        _poses.push_back(pos);
                }
        }

        inline Position findPosition(const Key& key) {
#ifdef PERFORMANCE_STATISTICS
                ++num_find_pos;
#endif
                Position pos = _hash(key) & (_ht_size - 1);

#ifdef PERFORMANCE_STATISTICS
                ++num_probes;
#endif
                if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                        return pos;
                }

                while (true) {
                        ++pos;
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        // the last element of _ht is always _empty_element
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                return pos;
                        }

                        ++pos;
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        // the last element of _ht is always _empty_element
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                return pos;
                        }

                        ++pos;
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        // the last element of _ht is always _empty_element
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                return pos;
                        }

                        ++pos;
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        // the last element of _ht is always _empty_element
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                return pos;
                        }

                        ++pos;
#ifdef PERFORMANCE_STATISTICS
                        ++num_probes;
#endif
                        // the last element of _ht is always _empty_element
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                return pos;
                        }
                }
        }

        Key _empty_element;
        uint64_t _max_size;
        uint64_t _ht_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
};

#ifdef PERFORMANCE_STATISTICS
template <typename Key, typename Value, typename Hash, bool TGrowable, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, SizeFactor>::num_access(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, SizeFactor>::num_contain(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, SizeFactor>::overall_max_size(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, SizeFactor>::num_probes(0);

template <typename Key, typename Value, typename Hash, bool TGrowable, size_t SizeFactor>
size_t HashMap<Key, Value, Hash, TGrowable, SizeFactor>::num_find_pos(0);
#endif

template <typename Key, typename Value, typename Hash = xxhash<Key>,
        bool TGrowable = false, bool Cache = true, size_t SizeFactor = 2>
class HashMapWithErase {
public:
        using Element = std::pair<Key, Value>;
        using key_type = Key;
        using mapped_type = Value;

private:
        using TSelf = HashMapWithErase<Key, Value, Hash, TGrowable, Cache, SizeFactor>;
        using Position = uint32_t;

public:
        using Iterator = HashTableWithEraseIterator<TSelf>;

        friend Iterator;

        explicit HashMapWithErase(const uint64_t max_size = 1) :
                _empty_element(std::numeric_limits<Key>::max()),
                _deleted_element(_empty_element - 1),
                _max_size(std::max(round_up_to_next_power_2(max_size), 16u)),
                _ht_size(max_size * SizeFactor),
                _size(0),
                _ht(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value())),
                _poses(),
                _hash(),
                _last_key(_empty_element),
                _last_position(0) {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(max_size);
        }

        static constexpr size_t size_factor = SizeFactor;

        HashMapWithErase(const TSelf&) = default;
        HashMapWithErase(TSelf&&) = default;

        TSelf& operator=(TSelf& other) = default;
        TSelf& operator=(TSelf&& other) = default;

        Iterator begin() const {
                return Iterator(*this, 0);
        }

        Iterator end() const {
                return Iterator(*this, _poses.size());
        }

        void reserve(const uint64_t max_size) {
                _max_size = round_up_to_next_power_2(max_size);
                _ht_size = _max_size * SizeFactor;
                _ht.resize(_ht_size + _max_size * 1.1, std::make_pair(_empty_element, Value()));
                _poses.reserve(_max_size);
        }

        inline size_t size() const {
                return _size;
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline void erase(const Key& key) {
                const Position pos = findPosition(key);
                if (_ht[pos].first != _empty_element && _ht[pos].first != _deleted_element) {
                        _ht[pos].first = _deleted_element;
                        --_size;
                }
        }

        inline Value& operator[] (const Key& key) {
                if (TGrowable && _poses.size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element || _ht[pos].first == _deleted_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = Value();
                        _poses.push_back(pos);
                        ++_size;
                }

                return _ht[pos].second;
        }

        inline bool contains(const Key& key) const {
                // Sometimes findPosition can return _delete_element
                // which means that in searched until the end of _ht and did not
                // find element but found _deleted_element. This means that
                // element is NOT in hash table
                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element || _ht[pos].first == _deleted_element) {
                        return false;
                }

                return true;
        }

        inline void insert(const Element& elem) {
                insertImpl(elem.first, elem.second);
        }

        inline void insert(const Key& key, const Value& value) {
                insertImpl(key, value);
        }

        inline void clear() {
                for (auto pos : _poses) {
                        _ht[pos].first = _empty_element;
                }
                _poses.clear();

                _last_key = _empty_element;
                _last_position = 0;
                _size = 0;
        }

        inline void swap(TSelf& hash_map) {
                std::swap(_ht_size, hash_map._ht_size);
                std::swap(_max_size, hash_map._max_size);
                _ht.swap(hash_map._ht);
                _poses.swap(hash_map._poses);
                std::swap(_last_key, hash_map._last_key);
                std::swap(_last_position, hash_map._last_position);
                std::swap(_empty_element, hash_map._empty_element);
                std::swap(_deleted_element, hash_map._deleted_element);
        }

private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos]);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key, const Value& value) {
                if (TGrowable && _poses.size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos].first == _empty_element || _ht[pos].first == _deleted_element) {
                        _ht[pos].first = key;
                        _ht[pos].second = value;
                        _poses.push_back(pos);
                        ++_size;
                }
        }

        inline Position findPosition(const Key& key) const {
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                const Position startPosition = _hash(key) & (_ht_size - 1);
                Position firstDeleted = std::numeric_limits<Position>::max();
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
                        if (_ht[pos].first == _empty_element || _ht[pos].first == key) {
                                if (firstDeleted != std::numeric_limits<Position>::max() && _ht[pos].first == key) {
                                        std::swap(_ht[firstDeleted], _ht[pos]);

                                        pos = firstDeleted;
                                }

                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }

                                return pos;
                        }
                        if (firstDeleted == std::numeric_limits<Position>::max() &&
                            _ht[pos].first == _deleted_element) {
                                firstDeleted = pos;
                        }
                }

                if (firstDeleted != std::numeric_limits<Position>::max()) {
                        if (Cache) {
                                _last_key = key;
                                _last_position = firstDeleted;
                        }
                        return firstDeleted;
                }

                std::cerr << "hash table overflowed" << std::endl;
                std::exit(-1);
        }

        key_type _empty_element;
        key_type _deleted_element;
        uint64_t _max_size;
        uint64_t _ht_size;
        uint64_t _size;
        mutable std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        mutable Key _last_key;
        mutable Position _last_position;
};

template <typename Key, typename Hash = xxhash<Key>, bool TGrowable = false,
        bool Cache = true, size_t SizeFactor = 2>
class HashSet {
public:
        using Element = Key;
        using hash_type = typename Hash::hash_type;
private:
        using TSelf = HashSet<Key, Hash, TGrowable, Cache, SizeFactor>;
        using Position = uint32_t;

        static constexpr hash_type max_hash_value = std::numeric_limits<hash_type>::max();
public:
        static constexpr size_t size_factor = SizeFactor;

        using Iterator = HashTableIterator<TSelf>;

        friend Iterator;

        explicit HashSet(const uint64_t max_size = 1) :
                _empty_element(std::numeric_limits<Key>::max()),
                _max_size(std::max(round_up_to_next_power_2(max_size), 16u)),
                _ht_size(_max_size * SizeFactor),
                _ht(_ht_size + _max_size * 1.1, _empty_element),
                _poses(),
                _hash(),
                _last_key(_empty_element),
                _last_position(0) {
                ALWAYS_ASSERT(is_power_2(_ht_size));
                _poses.reserve(max_size);
        }

        HashSet(const TSelf&) = default;
        HashSet(TSelf&&) = default;

        TSelf& operator= (const TSelf& other) {
                _ht_size = other._ht_size;
                _max_size = other._max_size;
                _ht = other._ht;
                _poses = other._poses;
                _hash = other._hash;
                _last_key = other._last_key;
                _last_position = other._last_position;
                return *this;
        }

        TSelf& operator= (TSelf&& other) {
                std::swap(_ht_size, other._ht_size);
                std::swap(_max_size, other._max_size);
                _ht.swap(other._ht);
                _poses.swap(other._poses);
                std::swap(_hash, other._hash);
                std::swap(_last_key, other._last_key);
                std::swap(_last_position, other._last_position);
                return *this;
        }

        void reserve(const uint32_t max_size) {
                _ht_size = max_size * SizeFactor;
                _max_size = max_size;
                _ht.resize(_ht_size + _max_size * 1.1, _empty_element);
                _poses.reserve(_max_size);
        }

        inline size_t size() const {
                return _poses.size();
        }

        Iterator begin() const {
                return Iterator(*this, 0);
        }

        Iterator end() const {
                return Iterator(*this, size());
        }

        inline bool empty() const {
                return size() == 0;
        }

        inline bool contains(const Key& key) {
                if (_ht[findPosition(key)] == _empty_element) {
                        return false;
                }
                return true;
        }

        inline void insert(const Key& key) {
                insertImpl(key);
        }

        inline void clear() {
                for (auto pos : _poses) {
                        _ht[pos] = _empty_element;
                }

                _poses.clear();

                _last_key = _empty_element;
                _last_position = 0;
        }

        inline void swap(TSelf& hash_set) {
                std::swap(_ht_size, hash_set._ht_size);
                std::swap(_max_size, hash_set._max_size);
                _ht.swap(hash_set._ht);
                _poses.swap(hash_set._poses);
                std::swap(_last_key, hash_set._last_key);
                std::swap(_last_position, hash_set._last_position);
        }

private:
        void resize() {
                TSelf new_hash_map(2 * _max_size);

                for (auto pos : _poses)
                        new_hash_map.insert(_ht[pos]);

                swap(new_hash_map);
        }

        inline void insertImpl(const Key& key) {
                if (TGrowable && size() == _max_size / 2)
                        resize();

                const Position pos = findPosition(key);
                if (_ht[pos] == _empty_element) {
                        _ht[pos] = key;
                        _poses.push_back(pos);
                }
        }

        inline Position findPosition(const Key& key) {
                if (Cache && key == _last_key) {
                        return _last_position;
                }

                //const Position startPosition = _hash(key) & (_ht_size - 1);
                const Position startPosition = _hash(key) % _ht_size;
                for (Position pos = startPosition; pos < _ht.size(); ++pos) {
                        if (_ht[pos] == _empty_element || _ht[pos] == key) {
                                if (Cache) {
                                        _last_key = key;
                                        _last_position = pos;
                                }
                                return pos;
                        }
                }

                throw std::runtime_error("Hash table overflowed");
        }

        const Key _empty_element;
        uint64_t _max_size;
        uint64_t _ht_size;
        std::vector<Element> _ht;
        std::vector<Position> _poses;
        Hash _hash;
        Key _last_key;
        Position _last_position;
};

//template <typename key_type, typename value_type>
//using hash_map = HashMap<key_type, value_type, simple_hash<key_type>, true, false>;

//template <typename key_type, typename value_type>
//using hash_map_with_erase = HashMapWithErase<key_type, value_type, simple_hash<key_type>, true>;
//
//template <typename key_type>
//using hash_set = HashSet<key_type, simple_hash<key_type>, true>;

//template <typename key_type, typename value_type>
//using hash_map = HashMap<key_type, value_type, xxhash<key_type>, true>;
//
//template <typename key_type, typename value_type>
//using hash_map_with_erase = HashMapWithErase<key_type, value_type, xxhash<key_type>, true>;
//
//template <typename key_type>
//using hash_set = HashSet<key_type, xxhash<key_type>, true>;

//template <typename key_type, typename value_type>
//using hash_map = HashMap<key_type, value_type, MurmurHash<key_type>, true>;

template <typename key_type, typename value_type>
using hash_map = HashMap<key_type, value_type, TabularHash<key_type, 5, 2, 10, true>, true>;

//template <typename key_type, typename value_type>
//using hash_map = HashMap<key_type, value_type, TabularHash<key_type, 0, 2, 10, true>, true>;

template <typename key_type, typename value_type>
using hash_map_with_erase = HashMapWithErase<key_type, value_type, MurmurHash<key_type>, true>;

template <typename key_type>
using hash_set = HashSet<key_type, MurmurHash<key_type>, true>;

}
