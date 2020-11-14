#pragma once

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "tools/macros_assertions.h"

#include <cstdint>

namespace parallel {

constexpr uint32_t g_cache_line_size = 64 * sizeof(char);

static size_t get_mem_for_thread(uint32_t proc_id, const PartitionConfig& config) {
        uint32_t num_threads = config.num_threads;
        uint32_t threads_per_socket = config.threads_per_socket;
        uint32_t l2_cache_size = config.l2_cache_size;
        uint32_t l3_cache_size = config.l3_cache_size;

        ALWAYS_ASSERT(proc_id < num_threads);
        uint32_t proc_in_full_socket = (num_threads / threads_per_socket) * threads_per_socket;
        uint32_t proc_per_socket;
        if (proc_id >= proc_in_full_socket) {
                proc_per_socket = num_threads - proc_in_full_socket;
        } else {
                proc_per_socket = threads_per_socket;
        }

        // two tables of size 1024 with uint32_t elements
        uint32_t tabular_hash_mem = 8 * 1024;
        ALWAYS_ASSERT(l3_cache_size > tabular_hash_mem * num_threads);
        return std::max(l3_cache_size / proc_per_socket, l2_cache_size - tabular_hash_mem);
}

static bool graph_is_large(graph_access& graph, const PartitionConfig& config) {
        return graph.number_of_edges() >= 1000000000;
}

template <typename T>
class alignas(g_cache_line_size) CacheAlignedData {
public:
        using Type = T;

        CacheAlignedData()
                :       m_elem()
        {
                static_assert(sizeof(Self) % g_cache_line_size == 0, "No cache line alignment");
        }

        template <typename... Args>
        CacheAlignedData(Args&&... args)
                :       m_elem(std::forward<Args>(args)...)
        {
                static_assert(sizeof(Self) % g_cache_line_size == 0, "No cache line alignment");
        }

        CacheAlignedData(CacheAlignedData& other) = default;
        CacheAlignedData(const CacheAlignedData& other) = default;
        CacheAlignedData(CacheAlignedData&& other) = default;

        inline T& get() {
                return m_elem;
        }

        inline const T& get() const {
                return m_elem;
        }

        inline bool operator! () const {
                return !m_elem;
        }

        inline CacheAlignedData& operator= (const T& elem) {
                m_elem = elem;
                return *this;
        }

        inline CacheAlignedData& operator= (const CacheAlignedData& other) {
                m_elem = other.m_elem;
                return *this;
        }

        inline CacheAlignedData& operator= (CacheAlignedData&& other) {
                m_elem = other.m_elem;
                return *this;
        }

        inline CacheAlignedData& operator+= (const T& elem) {
                m_elem += elem;
                return *this;
        }

        inline CacheAlignedData& operator-= (const T& elem) {
                m_elem -= elem;
                return *this;
        }

        inline CacheAlignedData& operator*= (const T& elem) {
                m_elem *= elem;
                return *this;
        }

        inline CacheAlignedData& operator/= (const T& elem) {
                m_elem /= elem;
                return *this;
        }

        inline T& operator++ () {
                ++m_elem;
                return *this;
        }

        inline T operator++ (int) {
                T tmp = m_elem;
                ++m_elem;
                return tmp;
        }
private:
        using Self = CacheAlignedData<T>;
        T m_elem;

};

}