#pragma once
#include "data_structure/graph_access.h"
#include "data_structure/parallel/graph_algorithm.h"
#include "data_structure/parallel/hash_table.h"
#include "data_structure/parallel/task_queue.h"
#include "data_structure/parallel/time.h"
#include "definitions.h"

#include <functional>
#include <vector>
#include <utility>

#include <tbb/concurrent_queue.h>

#include "data-structures/definitions.h"

namespace parallel {


class fast_boundary {
public:
        fast_boundary(graph_access& G, const PartitionConfig& config)
                :       m_G(G)
                ,       m_config(config)
                ,       m_blocks_info(m_G.get_partition_count())
        {}

        inline NodeWeight get_block_weight(PartitionID partition) {
                return m_blocks_info[partition].block_weight;
        }
        inline NodeID get_block_size(PartitionID partition) {
                return m_blocks_info[partition].block_size;
        }

        inline void set_block_weight(PartitionID partition, NodeWeight weight) {
                m_blocks_info[partition].block_weight = weight;
        }
        inline void set_block_size(PartitionID partition, NodeID size) {
                m_blocks_info[partition].block_size = size;
        }

        void balance_singletons() {
                for(size_t i = 0; i < m_singletons.size(); ++i) {
                        NodeWeight min = m_blocks_info[0].block_weight;
                        PartitionID p = 0;
                        for(size_t j = 0; j < m_blocks_info.size(); ++j) {
                                if (m_blocks_info[j].block_weight < min) {
                                        min = m_blocks_info[j].block_weight;
                                        p = j;
                                }
                        }

                        NodeID node = m_singletons[i];
                        if (m_blocks_info[p].block_weight + m_G.getNodeWeight(node) <= m_config.upper_bound_partition) {
                                m_blocks_info[m_G.getPartitionIndex(node)].block_weight -= m_G.getNodeWeight(node);
                                m_blocks_info[p].block_weight += m_G.getNodeWeight(node);
                                m_G.setPartitionIndex(node, p);
                        }
                }
        }

protected:
        using hash_map_with_erase_type = parallel::HashMapWithErase<NodeID, int32_t, parallel::xxhash<NodeID>, true, false>;

        struct block_data_type {
                NodeWeight block_weight = 0;
                NodeID block_size = 0;
        };

        graph_access& m_G;
        const PartitionConfig& m_config;
        std::vector<block_data_type> m_blocks_info;
        std::vector<NodeID> m_singletons;
};

class fast_sequential_boundary : public fast_boundary {
public:
        using iterator_type = hash_map_with_erase_type::Iterator;

        fast_sequential_boundary(graph_access& G, const PartitionConfig& config)
                :       fast_boundary(G, config)
                ,       m_boundary(1)
        {}

        iterator_type begin() const {
                return m_boundary.begin();
        }

        iterator_type end() const {
                return m_boundary.end();
        }

        void move(NodeID vertex, PartitionID from, PartitionID to) {
                ALWAYS_ASSERT(m_G.getPartitionIndex(vertex) == to);

                int32_t external_neighbors = 0;
                forall_out_edges(m_G, e, vertex) {
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID target_block = m_G.getPartitionIndex(target);

                        if (target_block == to) {
                                int32_t new_count = --m_boundary[target];
                                if (new_count == 0) {
                                        m_boundary.erase(target);
                                }
                        }

                        if (target_block == from) {
                                ++m_boundary[target];
                        }

                        if (target_block != to) {
                                ++external_neighbors;
                        }
                } endfor
                if (external_neighbors == 0) {
                        m_boundary.erase(vertex);
                } else {
                        m_boundary[vertex] = external_neighbors;
                }
        }

        void construct_boundary() {
                std::vector<std::pair<NodeID, int32_t>> preliminary_boundary;
                preliminary_boundary.reserve(1000);
                m_singletons.reserve(100);

                forall_nodes(m_G, n) {
                        int32_t num_external_neighbors = 0;

                        if (m_G.getNodeDegree(n) == 0) {
                                m_singletons.push_back(n);
                        }

                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        ++num_external_neighbors;
                                }
                        } endfor
                        if (num_external_neighbors > 0) {
                                preliminary_boundary.emplace_back(n, num_external_neighbors);
                        }

                        ++m_blocks_info[cur_block].block_size;
                        m_blocks_info[cur_block].block_weight += m_G.getNodeWeight(n);
                } endfor

                m_boundary.reserve(2 * preliminary_boundary.size());
                for (const auto& elem : preliminary_boundary) {
                        m_boundary[elem.first] = elem.second;
                }
        }

        void check_boundary() {
                std::vector<NodeID> this_boundary;
                this_boundary.reserve(m_boundary.size());

                for (const auto& elem : m_boundary) {
                        this_boundary.push_back(elem.first);
                }

                std::sort(this_boundary.begin(), this_boundary.end());

                std::vector<NodeID> expected_boundary;
                expected_boundary.reserve(m_boundary.size());
                std::cout << "Boundary diff: " << std::endl;
                forall_nodes(m_G, n) {
                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        bool boundary = false;
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        boundary = true;
                                        break;
                                }
                        } endfor
                        if (boundary) {
                                expected_boundary.push_back(n);
                                if (!m_boundary.contains(n)) {
                                        std::cout << "(" << n << ", FN) ";
                                }
                        }
                        if (!boundary) {
                                if (m_boundary.contains(n)) {
                                        std::cout << "(" << n << ", FP) ";
                                }
                        }
                } endfor
                std::cout << std::endl;
                std::cout << "check size = " << m_boundary.size() << std::endl;
                if (expected_boundary != this_boundary) {
                        std::cout << "expected size = " << expected_boundary.size() << std::endl;
                        std::cout << "this size = " << this_boundary.size() << std::endl;
                        ALWAYS_ASSERT(expected_boundary == this_boundary);
                }
        }
private:
        hash_map_with_erase_type m_boundary;

};

class fast_parallel_boundary : public fast_boundary {
private:
public:
        fast_parallel_boundary(graph_access& G, const PartitionConfig& config)
                :       fast_boundary(G, config)
                ,       m_primary_hash((uint32_t) m_config.seed)
                ,       m_secondary_hash((uint32_t) m_config.seed + 1)
                ,       m_second_level_size(m_config.num_threads * m_config.num_threads)
        {}

        // called only by single threaded executions
        void move(NodeID vertex, PartitionID from, PartitionID to) {
                ALWAYS_ASSERT(m_G.getPartitionIndex(vertex) == to);

                int32_t external_neighbors = 0;
                forall_out_edges(m_G, e, vertex) {
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID target_block = m_G.getPartitionIndex(target);

                        if (target_block == to) {
                                auto& target_boundary_ht = m_boundaries_per_thread[num_hash_table(target)].get();
                                int32_t new_count = --target_boundary_ht[target];
                                ALWAYS_ASSERT(new_count >= 0);
                                if (new_count == 0) {
                                        target_boundary_ht.erase(target);
                                }
                        }

                        if (target_block == from) {
                                auto& target_boundary_ht = m_boundaries_per_thread[num_hash_table(target)].get();
                                ++target_boundary_ht[target];
                        }

                        if (target_block != to) {
                                ++external_neighbors;
                        }
                } endfor

                auto& boundary_ht = m_boundaries_per_thread[num_hash_table(vertex)].get();
                if (external_neighbors == 0) {
                        boundary_ht.erase(vertex);
                } else {
                        boundary_ht[vertex] = external_neighbors;
                }
        }

        void construct_boundary() {
                container_collection_type containers(m_config.num_threads, Cvector<thread_container_type>(m_config.num_threads));

                if (m_config.num_threads == 1) {
                        distribute_boundary_vertices(m_G, containers, [this](uint32_t vertex, uint32_t) {
                                return external_neighbors(vertex);
                        });
                } else {
                        if (m_config.use_numa_aware_graph && parallel::graph_is_large(m_G, m_config)) {
                                parallel::numa_aware_graph na_graph(m_G, m_config.num_threads, m_config.threads_per_socket);
                                na_graph.construct();
                                distribute_boundary_vertices(
                                        na_graph.get_handle(0),
                                        containers,
                                        [this, &na_graph](uint32_t vertex, uint32_t thread_id) {
                                                return is_boundary(na_graph.get_handle(thread_id), vertex);
                                        }
                                );
                        } else {
                                distribute_boundary_vertices(m_G, containers, [this](uint32_t vertex, uint32_t) {
                                        return is_boundary(m_G, vertex);
                                });
                        }
                }

                // ------------------
                CLOCK_START;
                m_boundaries_per_thread.resize(m_config.num_threads);

                auto task_copy_to_hash_tables = [this, &containers] (uint32_t thread_id) {
                        NodeID total_size = 0;
                        auto& container = containers[thread_id];
                        for (const auto& sub_container : container) {
                                total_size += sub_container.get().size();
                        }

                        m_boundaries_per_thread[thread_id].get().reserve(std::max(total_size, 16u));

                        for (const auto& sub_container : container) {
                                for (const auto& elem : sub_container.get()) {
                                        m_boundaries_per_thread[thread_id].get().insert(elem.first, elem.second);
                                }
                        }
                        container.clear();
                        return total_size;
                };

                std::vector<std::future<NodeID>> futures_other;
                futures_other.reserve(parallel::g_thread_pool.NumThreads());
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures_other.push_back(parallel::g_thread_pool.Submit(i, task_copy_to_hash_tables, i + 1));
                }
                std::cout << task_copy_to_hash_tables(0) << " ";

                std::for_each(futures_other.begin(), futures_other.end(), [&](auto& future){
                        std::cout << future.get() << " ";
                });
                std::cout << std::endl;
                CLOCK_END("Copy to hash table");
        }

        const hash_map_with_erase_type& operator[] (uint32_t thread_id) const {
                return m_boundaries_per_thread[thread_id].get();
        }

        hash_map_with_erase_type& operator[] (uint32_t thread_id) {
                return m_boundaries_per_thread[thread_id].get();
        }

        void check_boundary() {
                std::vector<NodeID> this_boundary;

                size_t size = 0;
                for (uint32_t thread_id = 0; thread_id < m_config.num_threads; ++thread_id) {
                        size += m_boundaries_per_thread[thread_id].get().size();
                }
                this_boundary.reserve(size);

                for (uint32_t thread_id = 0; thread_id < m_config.num_threads; ++thread_id) {
                        for (const auto& elem : m_boundaries_per_thread[thread_id].get()) {
                                this_boundary.push_back(elem.first);
                        }
                }

                std::sort(this_boundary.begin(), this_boundary.end());

                std::vector<NodeID> expected_boundary;
                expected_boundary.reserve(size);
                std::cout << "Boundary diff: " << std::endl;
                forall_nodes(m_G, n) {
                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        bool boundary = false;
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        boundary = true;
                                        break;
                                }
                        } endfor

                        if (boundary) {
                                expected_boundary.push_back(n);
//                                if (!m_boundaries_per_thread[num_hash_table(n)].get().contains(n)) {
//                                     std::cout << "! " << n << " ";
//                                }
                        } else {
//                                if (m_boundaries_per_thread[num_hash_table(n)].get().contains(n)) {
//                                        std::cout << n << " ";
//                                }
                        }
                } endfor
                std::cout << std::endl;
                std::cout << "check size = " << size << std::endl;
                if (expected_boundary == this_boundary) {
                        std::cout << "expected size = " << expected_boundary.size() << std::endl;
                        std::cout << "this size = " << this_boundary.size() << std::endl;
                        ALWAYS_ASSERT(expected_boundary == this_boundary);
                }
        }

        ~fast_parallel_boundary() {
                if (m_boundaries_per_thread.empty()) {
                        return;
                }
                auto task = [this](uint32_t thread_id) {
                        m_boundaries_per_thread[thread_id].get().clear();
                        // to resize hash tables in parallel
                        m_boundaries_per_thread[thread_id].get().reserve(1);
                };
                std::vector<std::future<void>> futures;
                futures.reserve(parallel::g_thread_pool.NumThreads());
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, task, i + 1));
                }
                task(0);

                std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                        future.get();
                });
        }

        void begin_movements() {
                if (m_containers.size() != m_config.num_threads) {
                        m_containers.resize(m_config.num_threads, Cvector<thread_container<NodeID>>(m_second_level_size));
                }
        }

        void add_vertex_to_check(NodeID vertex) {
                uint32_t thread_id = num_hash_table(vertex);
                uint32_t bucket_id = (uint32_t) m_secondary_hash(vertex) % m_second_level_size;
                m_containers[thread_id][bucket_id].get().concurrent_emplace_back(vertex);
        }

        void finish_movements() {
                finish_movements_impl([&](uint32_t thread_id) -> hash_map_with_erase_type& {
                        return m_boundaries_per_thread[thread_id].get();
                });
        }
private:
        using thread_container_type = thread_container<std::pair<NodeID, int32_t>>;
        using container_collection_type = std::vector<Cvector<thread_container_type>>;


        Cvector<hash_map_with_erase_type> m_boundaries_per_thread;

        template <typename TGraph, typename TFunctor>
        void distribute_boundary_vertices(TGraph& graph, container_collection_type& containers, TFunctor&& bnd_func) {
                CLOCK_START;

                NodeID block_size = std::max<NodeID>(sqrt(graph.number_of_nodes()), 1000);

                std::atomic<NodeID> offset(0);
                auto task_distribute = [this, &containers, &offset, block_size, &bnd_func, &graph] (uint32_t thread_id) {
                        std::vector<block_data_type> blocks_info(m_G.get_partition_count());
                        NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                        while (begin < graph.number_of_nodes()) {
                                NodeID end = std::min(begin + block_size, graph.number_of_nodes());

                                for (NodeID node = begin; node != end; ++node) {
                                        PartitionID cur_part = graph.getPartitionIndex(node);
                                        int32_t num_external_neighbors = bnd_func(node, thread_id);

                                        ++blocks_info[cur_part].block_size;
                                        blocks_info[cur_part].block_weight += graph.getNodeWeight(node);
                                        if (num_external_neighbors > 0) {
                                                auto& container = containers[num_hash_table(node)][m_secondary_hash(node) % m_config.num_threads].get();
                                                container.concurrent_emplace_back(node, num_external_neighbors);
                                        }
                                }
                                begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                        }
                        return blocks_info;
                };

                m_blocks_info.clear();
                submit_for_all(task_distribute, [](auto& result, auto&& new_value) {
                        if (result.empty()) {
                                result = std::move(new_value);
                        } else {
                                for (size_t i = 0; i < result.size(); ++i) {
                                        result[i].block_size += new_value[i].block_size;
                                        result[i].block_weight += new_value[i].block_weight;
                                }
                        }
                }, m_blocks_info);
                CLOCK_END("Distribute boundary vertices");
        }

protected:
        std::vector<parallel::Cvector<parallel::thread_container<NodeID>>> m_containers;
        parallel::xxhash<NodeID> m_primary_hash;
        parallel::xxhash<NodeID> m_secondary_hash;
        uint32_t m_second_level_size;

        template <typename TFunctor>
        void finish_movements_impl(TFunctor&& get_ht) {
                auto task = [&, this](uint32_t thread_id) {
                        auto& container = m_containers[thread_id];
                        auto& ht = get_ht(thread_id);
                        parallel::HashSet<NodeID, parallel::xxhash<NodeID>, true> processed_vertices(1024u);
                        for (auto& sub_container : container) {
                                for (const auto& vertex : sub_container.get()) {
                                        if (processed_vertices.contains(vertex)) {
                                                continue;
                                        }
                                        processed_vertices.insert(vertex);
                                        int32_t num_external_neighbors = is_boundary(m_G, vertex);
                                        if (num_external_neighbors > 0) {
                                                // we do not care if num_external_neighbors
                                                // correct or not. It must be only greater than 0
                                                ht.insert(vertex, num_external_neighbors);
                                        } else {
                                                ht.erase(vertex);
                                        }
                                }
                                sub_container.get().clear();
                        }
                };
                parallel::submit_for_all(task);
        }

        inline NodeID num_hash_table(NodeID vertex) const {
                return (NodeID) m_primary_hash(vertex) % m_config.num_threads;
        }

        inline int32_t external_neighbors(NodeID vertex) const {
                int32_t num_external_neighbors = 0;
                PartitionID cur_part = m_G.getPartitionIndex(vertex);
                forall_out_edges(m_G, e, vertex){
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID part = m_G.getPartitionIndex(target);
                        if (cur_part != part) {
                                ++num_external_neighbors;
                        }
                } endfor
                return num_external_neighbors;
        }

        template <typename TGraph>
        inline int32_t is_boundary(TGraph& graph, NodeID vertex) const {
                int32_t num_external_neighbors = 0;
                PartitionID cur_part = graph.getPartitionIndex(vertex);
                forall_out_edges(graph, e, vertex) {
                        NodeID target = graph.getEdgeTarget(e);
                        PartitionID part = graph.getPartitionIndex(target);

                        if (cur_part == part) {
                                continue;
                        } else {
                                ++num_external_neighbors;
                                break;
                        }
                } endfor
                return num_external_neighbors;
        }
};

class fast_parallel_boundary_exp : public fast_parallel_boundary {
private:
        using concurrent_ht_type = growt::uaGrow<parallel::MurmurHash<uint64_t>>;
        using concurrent_ht_handle_type = typename concurrent_ht_type::Handle;

public:
        fast_parallel_boundary_exp(graph_access& G, const PartitionConfig& config)
                :       fast_parallel_boundary(G, config)
                ,       m_boundary(std::min<size_t>(G.number_of_edges(), 1000000))
        {
                m_ht_handles.reserve(m_config.num_threads);
                for (size_t thread_id = 0; thread_id < m_config.num_threads; ++thread_id) {
                        m_ht_handles.emplace_back(m_boundary.getHandle());
                }
        }

        // called only by single threaded executions
        void move(NodeID vertex, PartitionID from, PartitionID to) {
                ALWAYS_ASSERT(m_G.getPartitionIndex(vertex) == to);

                auto& ht_handle = m_ht_handles[0].get();
                int32_t external_neighbors = 0;
                forall_out_edges(m_G, e, vertex) {
                        NodeID target = m_G.getEdgeTarget(e);
                        PartitionID target_block = m_G.getPartitionIndex(target);

                        if (target_block == to) {
                                auto it = ht_handle.find(target);
                                ALWAYS_ASSERT(it != ht_handle.end());

                                int32_t new_count = (*it).second;
                                if (--new_count == 0) {
                                        ht_handle.erase(target);
                                } else {
                                        ht_handle.update(target, [](size_t& lhs, const size_t& rhs) {
                                                return lhs = rhs;
                                        }, new_count);
                                }
                        }

                        if (target_block == from) {
                                ht_handle.insertOrUpdate(target, 1, [](size_t& lhs) {
                                        return lhs += 1;
                                });
                        }

                        if (target_block != to) {
                                ++external_neighbors;
                        }
                } endfor

                if (external_neighbors == 0) {
                        ht_handle.erase(vertex);
                } else {
                        ht_handle.insertOrUpdate(vertex, external_neighbors, [](size_t& lhs, const size_t& rhs) {
                                return lhs = rhs;
                        }, external_neighbors);
                }
        }

        void construct_boundary() {
                if (parallel::g_thread_pool.NumThreads() == 0) {
                        distribute_boundary_vertices(m_G, [this](uint32_t vertex, uint32_t) {
                                return external_neighbors(vertex);
                        });
                } else {
                        if (m_config.use_numa_aware_graph && parallel::graph_is_large(m_G, m_config)) {
                                parallel::numa_aware_graph na_graph(m_G, m_config.num_threads, m_config.threads_per_socket);
                                na_graph.construct();
                                distribute_boundary_vertices(na_graph.get_handle(0), [this, &na_graph](uint32_t vertex, uint32_t thread_id) {
                                        return is_boundary(na_graph.get_handle(thread_id), vertex);
                                });
                        } else {
                                distribute_boundary_vertices(m_G, [this](uint32_t vertex, uint32_t) {
                                        return is_boundary(m_G, vertex);
                                });
                        }
                }
                std::cout << "Boundary size\t" << m_ht_handles[0].get().element_count_approx() << std::endl;
        }

        const concurrent_ht_handle_type& operator[] (uint32_t thread_id) const {
                return m_ht_handles[thread_id].get();
        }

        concurrent_ht_handle_type& operator[] (uint32_t thread_id) {
                return m_ht_handles[thread_id].get();
        }

        void check_boundary() {
                std::vector<NodeID> this_boundary;

                auto& ht_handle = m_ht_handles[0].get();
                this_boundary.reserve(ht_handle.element_count_approx());

                for (auto it = ht_handle.begin(); it != ht_handle.end(); ++it) {
                        this_boundary.push_back((*it).first);
                }
                size_t size = this_boundary.size();

                std::sort(this_boundary.begin(), this_boundary.end());

                std::vector<NodeID> expected_boundary;
                expected_boundary.reserve(size);
                std::cout << "Boundary diff: " << std::endl;
                forall_nodes(m_G, n) {
                        PartitionID cur_block = m_G.getPartitionIndex(n);
                        bool boundary = false;
                        forall_out_edges(m_G, e, n) {
                                NodeID target = m_G.getEdgeTarget(e);
                                PartitionID target_block = m_G.getPartitionIndex(target);

                                if (cur_block != target_block) {
                                        boundary = true;
                                        break;
                                }
                        } endfor

                        if (boundary) {
                                expected_boundary.push_back(n);
//                                if (!m_boundaries_per_thread[num_hash_table(n)].get().contains(n)) {
//                                     std::cout << "! " << n << " ";
//                                }
                        } else {
//                                if (m_boundaries_per_thread[num_hash_table(n)].get().contains(n)) {
//                                        std::cout << n << " ";
//                                }
                        }
                } endfor
                std::cout << std::endl;
                std::cout << "check size = " << size << std::endl;
                if (expected_boundary == this_boundary) {
                        std::cout << "expected size = " << expected_boundary.size() << std::endl;
                        std::cout << "this size = " << this_boundary.size() << std::endl;
                        ALWAYS_ASSERT(expected_boundary == this_boundary);
                }
        }

        ~fast_parallel_boundary_exp() {
        }

        //
        void begin_movements() {
                if (m_containers.size() != m_config.num_threads) {
                        m_containers.resize(m_config.num_threads, Cvector<thread_container<NodeID>>(m_second_level_size));
                }
        }

        void add_vertex_to_check(NodeID vertex) {
                uint32_t thread_id = num_hash_table(vertex);
                uint32_t bucket_id = (uint32_t) m_secondary_hash(vertex) % m_second_level_size;
                m_containers[thread_id][bucket_id].get().concurrent_emplace_back(vertex);
        }

        void finish_movements() {
                finish_movements_impl([&](uint32_t thread_id) -> concurrent_ht_handle_type& {
                        return m_ht_handles[thread_id].get();
                });
        }
private:
        concurrent_ht_type m_boundary;
        Cvector<concurrent_ht_handle_type> m_ht_handles;

        template <typename TGraph, typename TFunctor>
        void distribute_boundary_vertices(TGraph& graph, TFunctor&& bnd_func) {
                CLOCK_START;
                NodeID block_size = std::max<NodeID>(sqrt(graph.number_of_nodes()), 1000);

                std::atomic<NodeID> offset(0);
                auto task_distribute = [this, &offset, block_size, &bnd_func, &graph] (uint32_t thread_id) {
                        std::vector<block_data_type> blocks_info(m_G.get_partition_count());
                        NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                        auto& ht_handle = m_ht_handles[thread_id].get();
                        while (begin < graph.number_of_nodes()) {
                                NodeID end = std::min(begin + block_size, graph.number_of_nodes());

                                for (NodeID node = begin; node != end; ++node) {
                                        PartitionID cur_part = graph.getPartitionIndex(node);
                                        int32_t num_external_neighbors = bnd_func(node, thread_id);

                                        ++blocks_info[cur_part].block_size;
                                        blocks_info[cur_part].block_weight += graph.getNodeWeight(node);
                                        if (num_external_neighbors > 0) {
                                                ht_handle.insert(node, num_external_neighbors);
                                        }
                                }
                                begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                        }
                        return blocks_info;
                };

                m_blocks_info.clear();
                submit_for_all(task_distribute, [](auto& result, auto&& new_value) {
                        if (result.empty()) {
                                result = std::move(new_value);
                        } else {
                                for (size_t i = 0; i < result.size(); ++i) {
                                        result[i].block_size += new_value[i].block_size;
                                        result[i].block_weight += new_value[i].block_weight;
                                }
                        }
                }, m_blocks_info);
                CLOCK_END("Distribute boundary vertices");
        }
};

//using boundary_type = parallel::fast_sequential_boundary;
//using boundary_type = parallel::fast_parallel_boundary;
using boundary_type = parallel::fast_parallel_boundary_exp;

}