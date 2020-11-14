/******************************************************************************
 * contraction.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "contraction.h"
#include "data_structure/parallel/time.h"
#include "data_structure/parallel/thread_pool.h"
#include "../uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"
#include "macros_assertions.h"

#include <ips4o/ips4o.hpp>

#include "data-structures/definitions.h"

contraction::contraction() {

}

contraction::~contraction() {

}

// for documentation see technical reports of christian schulz  
void contraction::contract(const PartitionConfig& partition_config,
                           graph_access& G,
                           graph_access& coarser,
                           const Matching& edge_matching,
                           const CoarseMapping& coarse_mapping,
                           const NodeID& no_of_coarse_vertices,
                           const NodePermutationMap& permutation) const {
        coarser.setUnitWeightEdges(false);
        if (partition_config.matching_type == CLUSTER_COARSENING) {
                if (!partition_config.fast_contract_clustering) {
                        contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping,
                                            no_of_coarse_vertices, permutation);
                        return;
                } else {
//                        //return fast_contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping, no_of_coarse_vertices, permutation);
                        parallel_fast_contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping,
                                                          no_of_coarse_vertices, permutation);
                        return;
                }
        } else if (partition_config.matching_type == MATCHING_PARALLEL_LOCAL_MAX) {
                parallel_contract_matching(partition_config, G, coarser, edge_matching, coarse_mapping,
                                           no_of_coarse_vertices, permutation);
                return;
        }

        if (partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        std::vector<NodeID> new_edge_targets(G.number_of_edges());
        forall_edges(G, e){
                                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
                        }endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        //we dont know the number of edges jet, so we use the old number for 
        //construction of the coarser graph and then resize the field according
        //to the number of edges we really got
        coarser.start_construction(no_of_coarse_vertices, G.number_of_edges());

        NodeID cur_no_vertices = 0;

        forall_nodes(G, n){
                                NodeID node = permutation[n];
                                //we look only at the coarser nodes
                                if (coarse_mapping[node] != cur_no_vertices)
                                        continue;

                                NodeID coarseNode = coarser.new_node();
                                coarser.setNodeWeight(coarseNode, G.getNodeWeight(node));

                                if (partition_config.combine) {
                                        coarser.setSecondPartitionIndex(coarseNode, G.getSecondPartitionIndex(node));
                                }

                                // do something with all outgoing edges (in auxillary graph)
                                forall_out_edges(G, e, node){
                                                        visit_edge(G, coarser, edge_positions, coarseNode, e,
                                                                   new_edge_targets);
                                                }endfor

                                //this node was really matched
                                NodeID matched_neighbor = edge_matching[node];
                                if (node != matched_neighbor) {
                                        //update weight of coarser node
                                        NodeWeight new_coarse_weight =
                                                G.getNodeWeight(node) + G.getNodeWeight(matched_neighbor);
                                        coarser.setNodeWeight(coarseNode, new_coarse_weight);

                                        forall_out_edges(G, e, matched_neighbor){
                                                                visit_edge(G, coarser, edge_positions, coarseNode, e,
                                                                           new_edge_targets);
                                                        }endfor
                                }
                                forall_out_edges(coarser, e, coarseNode){
                                                        edge_positions[coarser.getEdgeTarget(e)] = UNDEFINED_EDGE;
                                                }endfor

                                cur_no_vertices++;
                        }endfor

        ASSERT_RANGE_EQ(edge_positions, 0, edge_positions.size(), UNDEFINED_EDGE);
        ASSERT_EQ(no_of_coarse_vertices, cur_no_vertices);

        //this also resizes the edge fields ... 
        coarser.finish_construction();
}

void contraction::contract_clustering(const PartitionConfig& partition_config,
                                      graph_access& G,
                                      graph_access& coarser,
                                      const Matching& edge_matching,
                                      const CoarseMapping& coarse_mapping,
                                      const NodeID& no_of_coarse_vertices,
                                      const NodePermutationMap& permutation) const {

        if (partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        //save partition map -- important if the graph is allready partitioned
        std::vector<int> partition_map(G.number_of_nodes());
        int k = G.get_partition_count();
        forall_nodes(G, node){
                                partition_map[node] = G.getPartitionIndex(node);
                                G.setPartitionIndex(node, coarse_mapping[node]);
                        }endfor

        G.set_partition_count(no_of_coarse_vertices);

        complete_boundary bnd(&G);
        bnd.build();
        bnd.getUnderlyingQuotientGraph(coarser);

        G.set_partition_count(k);
        forall_nodes(G, node){
                                G.setPartitionIndex(node, partition_map[node]);
                                coarser.setPartitionIndex(coarse_mapping[node], G.getPartitionIndex(node));

                                if (partition_config.combine) {
                                        coarser.setSecondPartitionIndex(coarse_mapping[node],
                                                                        G.getSecondPartitionIndex(node));
                                }

                        }endfor

}

void contraction::parallel_fast_contract_clustering(const PartitionConfig& partition_config,
                                                    graph_access& G,
                                                    graph_access& coarser,
                                                    const Matching&,
                                                    const CoarseMapping& coarse_mapping,
                                                    const NodeID& no_of_coarse_vertices,
                                                    const NodePermutationMap&) const {
        if (partition_config.balls_and_bins_ht) {
                parallel_fast_contract_clustering_multiple_threads_balls_and_bins_ht(partition_config, G,
                                                                                     coarser, coarse_mapping,
                                                                                     no_of_coarse_vertices);
                return;
        }

        CLOCK_START;
        // build set of new edges
        double avg_degree = (G.number_of_edges() + 0.0) / G.number_of_nodes();
        size_t num_cut_edges = std::min<size_t>(avg_degree * no_of_coarse_vertices, G.number_of_edges() / 10);
        std::cout << "overall ht capacity\t" << num_cut_edges << std::endl;

        growt::uaGrow<parallel::MurmurHash<uint64_t>> new_edges(2 * num_cut_edges);
        CLOCK_END("Init hash tables");

        CLOCK_START_N;
        std::atomic<NodeID> offset(0);

        NodeID block_size = (NodeID) sqrt(G.number_of_nodes());
        block_size = std::max(block_size, 1000u);
        std::cout << "block_size\t" << block_size << std::endl;

        auto task = [&](uint32_t id) {
                auto handle = new_edges.getHandle();
                std::vector<NodeWeight> my_block_infos(no_of_coarse_vertices);
                NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                while (begin < G.number_of_nodes()) {
                        NodeID end = std::min(begin + block_size, G.number_of_nodes());

                        for (NodeID node = begin; node != end; ++node) {
                                PartitionID source_cluster = coarse_mapping[node];
                                my_block_infos[source_cluster] += G.getNodeWeight(node);

                                forall_out_edges(G, e, node) {
                                        NodeID targetID = G.getEdgeTarget(e);
                                        PartitionID target_cluster = coarse_mapping[targetID];
                                        bool is_cut_edge = source_cluster != target_cluster;

                                        if (is_cut_edge) {
                                                EdgeWeight edge_weight = G.getEdgeWeight(e);
                                                uint64_t key = get_uint64_from_pair_unsorted(source_cluster,
                                                                                             target_cluster);

                                                handle.insertOrUpdate(key, edge_weight,
                                                                      [](size_t& lhs, const size_t& rhs) {
                                                                              return lhs += rhs;
                                                                      },
                                                                      edge_weight);
                                        }
                                } endfor
                        }
                        begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }
                return my_block_infos;
        };

        std::vector<NodeWeight> block_infos;
        parallel::submit_for_all(task, [&](auto& block_infos, auto&& cur_block_infos) {
                if (block_infos.empty()) {
                        block_infos.swap(cur_block_infos);
                } else {
                        for (size_t i = 0; i < no_of_coarse_vertices; ++i) {
                                block_infos[i] += cur_block_infos[i];
                        }
                }
        }, block_infos);
        CLOCK_END("Construct hash table and aux data");

        CLOCK_START_N;
        std::vector<parallel::AtomicWrapper<EdgeID>> offsets(no_of_coarse_vertices);
        offset.store(0, std::memory_order_relaxed);
        auto task1 = [&](uint32_t thread_id) {
                auto handle = new_edges.getHandle();
                size_t size = handle.capacity();
                const EdgeID block_size = std::max<EdgeID>(sqrt(size), 1000);
                EdgeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                EdgeID num_edges = 0;

                while (begin < size) {
                        auto it = handle.range(begin, begin + block_size);
                        for (; it != handle.range_end(); ++it) {
                                std::pair<NodeID, NodeID> edge = get_pair_from_uint64((*it).first);
                                ++num_edges;
                                offsets[edge.first].fetch_add(1, std::memory_order_relaxed);
                        }
                        begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }

                return num_edges;
        };


        EdgeID num_edges = parallel::submit_for_all(task1,
                                                    [](EdgeID num_edges, EdgeID cur_num_edges) {
                                                            return num_edges + cur_num_edges;
                                                    }, EdgeID(0));
        std::cout << "num edges\t" << num_edges << std::endl;
        CLOCK_END("Calculate offsets");

        CLOCK_START_N;
        NodeID last_node_degree = offsets.back();
        {
                CLOCK_START;
                parallel::partial_sum_open_interval(offsets.begin(), offsets.end(), offsets.begin(),
                                                    partition_config.num_threads);
                CLOCK_END("prefix sum");
        }
        std::vector<Node> nodes(no_of_coarse_vertices + 1);
        parallel::parallel_for_index(NodeID(0), no_of_coarse_vertices, [&](NodeID node) {
                nodes[node].firstEdge = offsets[node];
                nodes[node].weight = block_infos[node];
        });
        nodes.back().firstEdge = num_edges;
        ALWAYS_ASSERT(num_edges == offsets.back() + last_node_degree);
        CLOCK_END("Calculate prefix sum");

        CLOCK_START_N;
        std::vector<Edge> edges(num_edges);
        offset.store(0, std::memory_order_relaxed);
        auto task2 = [&](uint32_t thread_id) {
                auto handle = new_edges.getHandle();
                size_t size = handle.capacity();
                const EdgeID block_size = std::max<EdgeID>(sqrt(size), 1000);
                EdgeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);

                while (begin < size) {
                        auto it = handle.range(begin, begin + block_size);
                        for (; it != handle.range_end(); ++it) {
                                std::pair<NodeID, NodeID> edge = get_pair_from_uint64((*it).first);
                                auto edge_weight = (*it).second;
                                EdgeID offset = offsets[edge.first].fetch_add(1, std::memory_order_relaxed);
                                edges[offset].target = edge.second;
                                edges[offset].weight = edge_weight;
                        }
                        begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }
        };

        parallel::submit_for_all(task2);

        coarser.start_construction(nodes, edges);
        ALWAYS_ASSERT(!partition_config.graph_allready_partitioned);

        CLOCK_END("Calculate edges array");
}

void contraction::parallel_fast_contract_clustering_multiple_threads_balls_and_bins_ht(
        const PartitionConfig& partition_config,
        graph_access& G,
        graph_access& coarser,
        const CoarseMapping& coarse_mapping,
        const NodeID& no_of_coarse_vertices) const {
        CLOCK_START;
        const uint32_t num_threads = partition_config.num_threads;

        // build set of new edges
        double avg_degree = (G.number_of_edges() + 0.0) / G.number_of_nodes();
        size_t num_cut_edges = std::min<size_t>(avg_degree * no_of_coarse_vertices, G.number_of_edges() / 10);
        std::cout << "overall ht capacity\t" << num_cut_edges << std::endl;
        std::cout << "ht capacity\t" << num_cut_edges / num_threads << std::endl;

        using concurrent_ht_type = growt::uaGrow<parallel::MurmurHash<uint64_t>>;
        parallel::ParallelVector<concurrent_ht_type> new_edges(num_threads);

        parallel::submit_for_all([&](uint32_t thread_id) {
                new(&new_edges[thread_id]) concurrent_ht_type(2 * num_cut_edges / num_threads);
        });
        CLOCK_END("Init hash tables");

        CLOCK_START_N;
        std::atomic<NodeID> offset(0);

        NodeID block_size = (NodeID) sqrt(G.number_of_nodes());
        block_size = std::max(block_size, 1000u);
        std::cout << "block_size\t" << block_size << std::endl;

        auto task_with_buffers = [&](uint32_t) {
                std::vector<concurrent_ht_type::Handle> handles;
                handles.reserve(num_threads);
                for (size_t i = 0; i < num_threads; ++i) {
                        handles.push_back(new_edges[i].getHandle());
                }
                std::vector<NodeWeight> my_block_infos(no_of_coarse_vertices);
                std::vector<std::vector<std::pair<uint64_t, EdgeWeight>>> buffers(num_threads);
                const NodeID max_buffer_size = 10000;
                for (auto& buffer : buffers) {
                        buffer.reserve(max_buffer_size);
                }

                auto empty_buffer_task = [&](auto& handle, auto& buffer) {
                        for (const auto& elem : buffer) {
                                uint64_t key = elem.first;
                                EdgeWeight edge_weight = elem.second;
                                handle.insertOrUpdate(key, edge_weight,
                                                      [](size_t& lhs, const size_t& rhs) {
                                                              return lhs += rhs;
                                                      },
                                                      edge_weight);
                        }
                        buffer.clear();
                };

                NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                while (begin < G.number_of_nodes()) {
                        NodeID end = std::min(begin + block_size, G.number_of_nodes());

                        for (NodeID node = begin; node != end; ++node) {
                                const PartitionID source_cluster = coarse_mapping[node];
                                my_block_infos[source_cluster] += G.getNodeWeight(node);
                                uint32_t ht_num = num_threads;

                                forall_out_edges(G, e, node){
                                                        NodeID targetID = G.getEdgeTarget(e);
                                                        PartitionID target_cluster = coarse_mapping[targetID];
                                                        bool is_cut_edge = source_cluster != target_cluster;

                                                        if (is_cut_edge) {
                                                                if (ht_num == num_threads) {
                                                                        ht_num = source_cluster % ht_num;
                                                                }
                                                                EdgeWeight edge_weight = G.getEdgeWeight(e);
                                                                uint64_t key = get_uint64_from_pair_unsorted(
                                                                        source_cluster, target_cluster);
                                                                buffers[ht_num].emplace_back(key, edge_weight);
                                                        }
                                                }endfor

                                if (ht_num != num_threads && buffers[ht_num].size() >= max_buffer_size) {
                                        empty_buffer_task(handles[ht_num], buffers[ht_num]);
                                }
                        }
                        begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }

                for (size_t i = 0; i < num_threads; ++i) {
                        empty_buffer_task(handles[i], buffers[i]);
                }

                return my_block_infos;
        };

        auto task_without_buffers = [&](uint32_t) {
                std::vector<concurrent_ht_type::Handle> handles;
                handles.reserve(num_threads);
                for (size_t i = 0; i < num_threads; ++i) {
                        handles.push_back(new_edges[i].getHandle());
                }
                std::vector<NodeWeight> my_block_infos(no_of_coarse_vertices);

                NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                while (begin < G.number_of_nodes()) {
                        NodeID end = std::min(begin + block_size, G.number_of_nodes());

                        for (NodeID node = begin; node != end; ++node) {
                                const PartitionID source_cluster = coarse_mapping[node];
                                my_block_infos[source_cluster] += G.getNodeWeight(node);
                                concurrent_ht_type::Handle* handle_ptr = nullptr;

                                forall_out_edges(G, e, node){
                                        NodeID targetID = G.getEdgeTarget(e);
                                        PartitionID target_cluster = coarse_mapping[targetID];
                                        bool is_cut_edge = source_cluster != target_cluster;

                                        if (is_cut_edge) {
                                                EdgeWeight edge_weight = G.getEdgeWeight(e);
                                                uint64_t key = get_uint64_from_pair_unsorted(source_cluster,
                                                                                             target_cluster);

                                                if (handle_ptr == nullptr) {
                                                        handle_ptr = &handles[source_cluster % num_threads];
                                                }

                                                handle_ptr->insertOrUpdate(key, edge_weight,
                                                                           [](size_t& lhs, const size_t& rhs) {
                                                                                   return lhs += rhs;
                                                                           }, edge_weight);
                                        }
                                } endfor
                        }
                        begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }
                return my_block_infos;
        };

        auto task = task_without_buffers;

        std::vector<NodeWeight> block_infos;
        block_infos.reserve(no_of_coarse_vertices);
        parallel::submit_for_all(task, [&](auto& block_infos, auto&& cur_block_infos) {
                if (block_infos.empty()) {
                        block_infos.swap(cur_block_infos);
                } else {
                        for (size_t i = 0; i < no_of_coarse_vertices; ++i) {
                                block_infos[i] += cur_block_infos[i];
                        }
                }
        }, block_infos);

        CLOCK_END("Construct hash table and aux data");

        CLOCK_START_N;
        std::vector<EdgeID> offsets(no_of_coarse_vertices);
        auto task1 = [&](uint32_t thread_id) {
                auto handle = new_edges[thread_id].getHandle();
                EdgeID num_edges = 0;
                for (auto it = handle.begin(); it != handle.end(); ++it) {
                        std::pair<NodeID, NodeID> edge = get_pair_from_uint64((*it).first);
                        ++num_edges;
                        ++offsets[edge.first];
                }
                return num_edges;
        };


        EdgeID num_edges = parallel::submit_for_all(task1,
                                                    [](EdgeID num_edges, EdgeID cur_num_edges) {
                                                            return num_edges + cur_num_edges;
                                                    }, EdgeID(0));
        std::cout << "num edges\t" << num_edges << std::endl;
        CLOCK_END("Calculate offsets");

        CLOCK_START_N;
        NodeID last_node_degree = offsets.back();
        {
                CLOCK_START;
                parallel::partial_sum_open_interval(offsets.begin(), offsets.end(), offsets.begin(),
                                                    partition_config.num_threads);
                CLOCK_END("prefix sum");
        }
        std::vector<Node> nodes(no_of_coarse_vertices + 1);
        parallel::parallel_for_index(NodeID(0), no_of_coarse_vertices, [&](NodeID node) {
                nodes[node].firstEdge = offsets[node];
                nodes[node].weight = block_infos[node];
        });
        nodes.back().firstEdge = num_edges;
        ALWAYS_ASSERT(num_edges == offsets.back() + last_node_degree);
        CLOCK_END("Calculate prefix sum");

        CLOCK_START_N;
        std::vector<Edge> edges(num_edges);
        auto task2 = [&](uint32_t thread_id) {
                auto handle = new_edges[thread_id].getHandle();
                for (auto it = handle.begin(); it != handle.end(); ++it) {
                        std::pair < NodeID, NodeID > edge = get_pair_from_uint64((*it).first);
                        auto edge_weight = (*it).second;
                        edges[offsets[edge.first]].target = edge.second;
                        edges[offsets[edge.first]].weight = edge_weight;
                        ++offsets[edge.first];
                }
        };

        parallel::submit_for_all(task2);

        coarser.start_construction(nodes, edges);
        ALWAYS_ASSERT(!partition_config.graph_allready_partitioned);

        CLOCK_END("Calculate edges array");

        CLOCK_START_N;
        auto task_clean_ht = [&](uint32_t thread_id) {
                new_edges[thread_id].~concurrent_ht_type();
        };
        parallel::submit_for_all(task_clean_ht);
        CLOCK_END("Clean hash tables");
}

void contraction::parallel_contract_matching(const PartitionConfig& partition_config,
                                             graph_access& G,
                                             graph_access& coarser,
                                             const Matching& edge_matching,
                                             const CoarseMapping& coarse_mapping,
                                             const NodeID& no_of_coarse_vertices,
                                             const NodePermutationMap&) const {
        NodeID node_block_size = (NodeID) sqrt(G.number_of_nodes());
        node_block_size = std::max(node_block_size, 1000u);

        CLOCK_START;
        std::vector<EdgeID> offsets(no_of_coarse_vertices);
        std::vector<NodeWeight> coarse_node_weights(no_of_coarse_vertices);
        std::atomic<NodeID> offset(0);
        auto task1 = [&](uint32_t thread_id) {
                parallel::hash_set<NodeID> common_neighbors(512);
                EdgeID num_edges = 0;
                while (true) {
                        NodeID begin = offset.fetch_add(node_block_size, std::memory_order_relaxed);
                        NodeID end = begin + node_block_size;
                        end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                        if (begin >= G.number_of_nodes()) {
                                break;
                        }

                        for (NodeID node = begin; node != end; ++node) {
                                NodeID matched = edge_matching[node];
                                EdgeID degree = 0;
                                if (node < matched) {
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (target != matched) {
                                                        common_neighbors.insert(coarse_mapping[target]);
                                                }
                                        } endfor

                                        forall_out_edges(G, e, matched) {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (target != node) {
                                                        common_neighbors.insert(coarse_mapping[target]);
                                                }
                                        } endfor

                                        coarse_node_weights[coarse_mapping[node]] = G.getNodeWeight(node) + G.getNodeWeight(matched);
                                } else if (node == matched) {
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                common_neighbors.insert(coarse_mapping[target]);
                                        } endfor
                                        coarse_node_weights[coarse_mapping[node]] = G.getNodeWeight(node);
                                } else {
                                        continue;
                                }
                                degree = common_neighbors.size();
                                common_neighbors.clear();
                                offsets[coarse_mapping[node]] = degree;
                                num_edges += degree;
                        }
                }
                return num_edges;
        };

        EdgeID num_edges = 0;
        {
                CLOCK_START;
                num_edges = parallel::submit_for_all(task1,
                                                     [](EdgeID num_edges, EdgeID cur_num_edges) {
                                                             return num_edges + cur_num_edges;
                                                     }, EdgeID(0));
                CLOCK_END("par");
        }
        CLOCK_END("Calculate offsets");


        CLOCK_START_N;
        NodeID last_node_degree = offsets.back();
        {
                CLOCK_START;
                parallel::partial_sum_open_interval(offsets.begin(), offsets.end(), offsets.begin(),
                                                    partition_config.num_threads);
                CLOCK_END("prefix sum");
        }
        std::vector<Node> nodes(no_of_coarse_vertices + 1);
        parallel::parallel_for_index(NodeID(0), no_of_coarse_vertices, [&](NodeID node) {
                nodes[node].firstEdge = offsets[node];
                nodes[node].weight = coarse_node_weights[node];
        });
        nodes.back().firstEdge = num_edges;
        ALWAYS_ASSERT(num_edges == offsets.back() + last_node_degree);
        CLOCK_END("Calculate prefix sum");

//        CLOCK_START_N;
//        std::vector<Node> nodes(no_of_coarse_vertices + 1);
//        {
//                CLOCK_START;
//                EdgeID cur_prefix = 0;
//                for (size_t i = 0; i < no_of_coarse_vertices; ++i) {
//                        EdgeID cur_offset = offsets[i];
//
//                        offsets[i] = cur_prefix;
//                        nodes[i].firstEdge = cur_prefix;
//                        nodes[i].weight = coarse_node_weights[i];
//
//                        cur_prefix += cur_offset;
//                }
//                nodes.back().firstEdge = cur_prefix;
//                CLOCK_END("prefix sum");
//        }
//        CLOCK_END("Calculate prefix sum");

        CLOCK_START_N;
        std::vector<Edge> edges(num_edges);
        CLOCK_END("Init edge array");

        offset.store(0, std::memory_order_relaxed);
        auto task2 = [&](uint32_t thread_id) {
                parallel::hash_map<NodeID, EdgeWeight> common_neighbors(512);
                while (true) {
                        NodeID begin = offset.fetch_add(node_block_size, std::memory_order_relaxed);
                        NodeID end = begin + node_block_size;
                        end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                        if (begin >= G.number_of_nodes()) {
                                break;
                        }

                        for (NodeID node = begin; node != end; ++node) {
                                NodeID matched = edge_matching[node];
                                if (node < matched) {
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (target != matched) {
                                                        common_neighbors[coarse_mapping[target]] += G.getEdgeWeight(e);
                                                }
                                        }endfor

                                        forall_out_edges(G, e, matched) {
                                                NodeID target = G.getEdgeTarget(e);
                                                if (target != node) {
                                                        common_neighbors[coarse_mapping[target]] += G.getEdgeWeight(e);
                                                }
                                        }endfor
                                } else if (node == matched) {
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                common_neighbors[coarse_mapping[target]] += G.getEdgeWeight(e);
                                        } endfor
                                } else {
                                        continue;
                                }

                                NodeID coarse_node = coarse_mapping[node];
                                for (const auto& record : common_neighbors) {
                                        edges[offsets[coarse_node]].target = record.first;
                                        edges[offsets[coarse_node]].weight = record.second;
                                        ++offsets[coarse_node];
                                }
                                common_neighbors.clear();
                        }
                }
        };

        CLOCK_START_N;
        parallel::submit_for_all(task2);
        CLOCK_END("Make edge array");

        CLOCK_START_N;
        coarser.start_construction(nodes, edges);
        CLOCK_END("Make graph");
}

void contraction::fast_contract_clustering(const PartitionConfig& partition_config,
                                           graph_access& G,
                                           graph_access& coarser,
                                           const Matching&,
                                           const CoarseMapping& coarse_mapping,
                                           const NodeID& no_of_coarse_vertices,
                                           const NodePermutationMap&) const {
        if (!partition_config.parallel_coarsening_lp && partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        CLOCK_START;
        std::vector<NodeWeight> block_infos(no_of_coarse_vertices);

        // build set of new edges
        double avg_degree = (G.number_of_edges() + 0.0) / G.number_of_nodes();
        size_t num_cut_edges = std::min<size_t>(avg_degree * no_of_coarse_vertices, G.number_of_edges() / 2);
        parallel::HashMap<uint64_t, EdgeWeight, parallel::MurmurHash<uint64_t>, true> new_edges(num_cut_edges);

        forall_nodes(G, n) {
                PartitionID source_cluster = coarse_mapping[n];
                block_infos[source_cluster] += G.getNodeWeight(n);

                forall_out_edges(G, e, n){
                        NodeID targetID = G.getEdgeTarget(e);
                        PartitionID target_cluster = coarse_mapping[targetID];
                        bool is_cut_edge = source_cluster != target_cluster;

                        if (is_cut_edge) {
                                new_edges[get_uint64_from_pair_sorted(source_cluster, target_cluster)] +=
                                        G.getEdgeWeight(e);
                        }
                } endfor
        } endfor
        CLOCK_END("Construct hash table and aux data");

        // construct graph
        CLOCK_START_N;
        std::vector<std::vector<std::pair<PartitionID, EdgeWeight>>> building_tool(no_of_coarse_vertices);
        for (auto& data : building_tool) {
                data.reserve(avg_degree);
        }

        for (const auto& edge_data : new_edges) {
                auto edge = get_pair_from_uint64(edge_data.first);
                auto edge_weight = edge_data.second / 2;
                building_tool[edge.first].emplace_back(edge.second, edge_weight);
                building_tool[edge.second].emplace_back(edge.first, edge_weight);
        }

        coarser.start_construction(building_tool.size(), 2 * new_edges.size());

        for (size_t p = 0; p < building_tool.size(); p++) {
                NodeID node = coarser.new_node();
                coarser.setNodeWeight(node, block_infos[p]);

                for (size_t j = 0; j < building_tool[p].size(); j++) {
                        EdgeID e = coarser.new_edge(node, building_tool[p][j].first);
                        coarser.setEdgeWeight(e, building_tool[p][j].second);
                }
        }
        coarser.finish_construction();

        if (!partition_config.parallel_coarsening_lp) {
                forall_nodes(G, node) {
                        coarser.setPartitionIndex(coarse_mapping[node], G.getPartitionIndex(node));

                        if (partition_config.combine) {
                                coarser.setSecondPartitionIndex(coarse_mapping[node], G.getSecondPartitionIndex(node));
                        }
                } endfor
        }
        CLOCK_END("Construct graph");
}

// for documentation see technical reports of christian schulz  
void contraction::contract_partitioned(const PartitionConfig& partition_config,
                                       graph_access& G,
                                       graph_access& coarser,
                                       const Matching& edge_matching,
                                       const CoarseMapping& coarse_mapping,
                                       const NodeID& no_of_coarse_vertices,
                                       const NodePermutationMap& permutation) const {

        if (partition_config.matching_type == CLUSTER_COARSENING) {
                return contract_clustering(partition_config, G, coarser, edge_matching, coarse_mapping,
                                           no_of_coarse_vertices, permutation);
        }


        std::vector<NodeID> new_edge_targets(G.number_of_edges());
        forall_edges(G, e){
                                new_edge_targets[e] = coarse_mapping[G.getEdgeTarget(e)];
                        }endfor

        std::vector<EdgeID> edge_positions(no_of_coarse_vertices, UNDEFINED_EDGE);

        //we dont know the number of edges jet, so we use the old number for 
        //construction of the coarser graph and then resize the field according
        //to the number of edges we really got
        coarser.set_partition_count(G.get_partition_count());
        coarser.start_construction(no_of_coarse_vertices, G.number_of_edges());

        if (partition_config.combine) {
                coarser.resizeSecondPartitionIndex(no_of_coarse_vertices);
        }

        NodeID cur_no_vertices = 0;

        PRINT(std::cout << "contracting a partitioned graph" << std::endl;)
        forall_nodes(G, n){
                                NodeID node = permutation[n];
                                //we look only at the coarser nodes
                                if (coarse_mapping[node] != cur_no_vertices)
                                        continue;

                                NodeID coarseNode = coarser.new_node();
                                coarser.setNodeWeight(coarseNode, G.getNodeWeight(node));
                                coarser.setPartitionIndex(coarseNode, G.getPartitionIndex(node));

                                if (partition_config.combine) {
                                        coarser.setSecondPartitionIndex(coarseNode, G.getSecondPartitionIndex(node));
                                }
                                // do something with all outgoing edges (in auxillary graph)
                                forall_out_edges(G, e, node){
                                                        visit_edge(G, coarser, edge_positions, coarseNode, e,
                                                                   new_edge_targets);
                                                }endfor

                                //this node was really matched
                                NodeID matched_neighbor = edge_matching[node];
                                if (node != matched_neighbor) {
                                        //update weight of coarser node
                                        NodeWeight new_coarse_weight =
                                                G.getNodeWeight(node) + G.getNodeWeight(matched_neighbor);
                                        coarser.setNodeWeight(coarseNode, new_coarse_weight);

                                        forall_out_edges(G, e, matched_neighbor){
                                                                visit_edge(G, coarser, edge_positions, coarseNode, e,
                                                                           new_edge_targets);
                                                        }endfor
                                }
                                forall_out_edges(coarser, e, coarseNode){
                                                        edge_positions[coarser.getEdgeTarget(e)] = UNDEFINED_EDGE;
                                                }endfor

                                cur_no_vertices++;
                        }endfor

        ASSERT_RANGE_EQ(edge_positions, 0, edge_positions.size(), UNDEFINED_EDGE);
        ASSERT_EQ(no_of_coarse_vertices, cur_no_vertices);

        //this also resizes the edge fields ... 
        coarser.finish_construction();
}
