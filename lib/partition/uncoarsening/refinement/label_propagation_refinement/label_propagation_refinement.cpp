/******************************************************************************
 * label_propagation_refinement.cpp 
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

#include "data_structure/parallel/thread_pool.h"
#include "data_structure/parallel/time.h"
#include "label_propagation_refinement.h"
#include "partition/coarsening/clustering/node_ordering.h"

using namespace parallel;

label_propagation_refinement::label_propagation_refinement() {
}

label_propagation_refinement::~label_propagation_refinement() {
}

EdgeWeight label_propagation_refinement::perform_refinement(PartitionConfig& config, graph_access& G, complete_boundary&) {
        std::vector<uint8_t> boundary_nodes;
        if (!config.parallel_lp) {
                CLOCK_START;
                auto res = sequential_label_propagation(config, G);
                CLOCK_END("Uncoarsening: Sequential lp");

                return res;
        } else {
                CLOCK_START;
                auto res = parallel_label_propagation(config, G);
                CLOCK_END("Uncoarsening: Parallel lp");

                return res;
        }
}

EdgeWeight label_propagation_refinement::perform_refinement(PartitionConfig& config, graph_access& G) {

        if (!config.parallel_lp) {
                CLOCK_START;
                auto res = sequential_label_propagation(config, G);
                CLOCK_END("Uncoarsening: Sequential lp");

                return res;
        } else {
                CLOCK_START;
                auto res = parallel_label_propagation(config, G);
                CLOCK_END("Uncoarsening: Parallel lp");

                return res;
        }
}

EdgeWeight label_propagation_refinement::sequential_label_propagation(PartitionConfig & partition_config,
                                                                      graph_access & G) {
        auto begin = std::chrono::high_resolution_clock::now();
        NodeWeight block_upperbound = partition_config.upper_bound_partition;

        // in this case the _matching paramter is not used
        // coarse_mappng stores cluster id and the mapping (it is identical)
        std::vector<PartitionID> hash_map(partition_config.k,0);
        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(partition_config.k, 0);

        auto b = std::chrono::high_resolution_clock::now();
        node_ordering n_ordering;
        n_ordering.order_nodes(partition_config, G, permutation);
        auto e = std::chrono::high_resolution_clock::now();
        std::cout << "Uncoarsening: Sequential init of permutations lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()
                  << std::endl;

        std::queue< NodeID > * Q             = new std::queue< NodeID >();
        std::queue< NodeID > * next_Q        = new std::queue< NodeID >();
        std::vector<bool> * Q_contained      = new std::vector<bool>(G.number_of_nodes(), false);
        std::vector<bool> * next_Q_contained = new std::vector<bool> (G.number_of_nodes(), false);
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
                Q->push(permutation[node]);
        } endfor
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Uncoarsening: Init sequential lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;

        begin = std::chrono::high_resolution_clock::now();

        int parts_size = partition_config.k * 2;
        std::vector<PartitionID> partIDs;
        partIDs.reserve(parts_size);
        unsigned int change_counter = 0;

        for( int j = 0; j < partition_config.label_iterations_refinement; j++) {
                while( !Q->empty() ) {
                        NodeID node = Q->front();
                        Q->pop();
                        (*Q_contained)[node] = false;

                        if (G.getNodeDegree(node) <= parts_size) {
                                //now move the node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        PartitionID part = G.getPartitionIndex(target);
                                        hash_map[part] += G.getEdgeWeight(e);
                                        partIDs.push_back(part);
                                } endfor
                        } else {
                                //now move the node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                                } endfor
                        }

                        //second sweep for finding max and resetting array
                        PartitionID max_block = G.getPartitionIndex(node);
                        PartitionID my_block  = G.getPartitionIndex(node);
                        PartitionID max_value = 0;

                        if (G.getNodeDegree(node) <= parts_size) {
                                for (auto cur_block : partIDs) {
                                        if (hash_map[cur_block] == 0) {
                                                continue;
                                        }

                                        PartitionID cur_value     = hash_map[cur_block];

                                        if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool()))
                                           && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || (cur_block == my_block && cluster_sizes[my_block] <= partition_config.upper_bound_partition)))
                                        {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                        }

                                        hash_map[cur_block] = 0;
                                }
                                partIDs.clear();
                        } else {
                                forall_out_edges(G, e, node) {
                                        NodeID target             = G.getEdgeTarget(e);
                                        PartitionID cur_block     = G.getPartitionIndex(target);
                                        if (hash_map[cur_block] == 0) {
                                                continue;
                                        }

                                        PartitionID cur_value     = hash_map[cur_block];

                                        if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool()))
                                           && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || (cur_block == my_block && cluster_sizes[my_block] <= partition_config.upper_bound_partition)))
                                        {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                        }

                                        hash_map[cur_block] = 0;
                                } endfor
                        }

                        cluster_sizes[G.getPartitionIndex(node)]  -= G.getNodeWeight(node);
                        cluster_sizes[max_block]         += G.getNodeWeight(node);
                        bool changed_label                = G.getPartitionIndex(node) != max_block;
                        change_counter                   += changed_label;
                        G.setPartitionIndex(node, max_block);
                        //std::cout <<  "maxblock " <<  max_block  << std::endl;

                        if(changed_label) {
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        if(!(*next_Q_contained)[target]) {
                                                next_Q->push(target);
                                                (*next_Q_contained)[target] = true;
                                        }
                                } endfor
                        }
                }

                std::swap( Q, next_Q);
                std::swap( Q_contained, next_Q_contained);

        }

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Uncoarsening: Main sequential lp:\t" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()
                  << std::endl;
        std::cout << "Uncoarsening: Improved:\t" << change_counter << std::endl;

        delete Q;
        delete next_Q;
        delete Q_contained;
        delete next_Q_contained;

        // in this case the _matching paramter is not used
        // coarse_mappng stores cluster id and the mapping (it is identical)
        //std::vector<PartitionID> hash_map(G.number_of_nodes(),0);
        //std::vector<NodeID> permutation(G.number_of_nodes());
        //std::vector<NodeWeight> cluster_sizes(partition_config.k,0);

        //forall_nodes(G, node) {
                //cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        //} endfor

        //random_functions::permutate_vector_fast(permutation, true);
        //NodeWeight block_upperbound = partition_config.upper_bound_partition;

        //for( int j = 0; j < partition_config.label_iterations; j++) {
                //forall_nodes(G, i) {
                        //NodeID node = permutation[i];
                        ////move the node to the cluster that is most common in the neighborhood

                        //forall_out_edges(G, e, node) {
                                //NodeID target = G.getEdgeTarget(e);
                                //hash_map[G.getPartitionIndex(target)]+=G.getEdgeWeight(e);
                        //} endfor

                        ////second sweep for finding max and resetting array
                        //PartitionID max_block = G.getPartitionIndex(node);
                        //PartitionID my_block  = G.getPartitionIndex(node);

                        //PartitionID max_value = 0;
                        //forall_out_edges(G, e, node) {
                                //NodeID target             = G.getEdgeTarget(e);
                                //PartitionID cur_block     = G.getPartitionIndex(target);
                                //PartitionID cur_value     = hash_map[cur_block];
                                //if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool()))
                                //&& (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || cur_block == my_block))
                                //{
                                        //max_value = cur_value;
                                        //max_block = cur_block;
                                //}

                                //hash_map[cur_block] = 0;
                        //} endfor
                        //cluster_sizes[G.getPartitionIndex(node)] -= G.getNodeWeight(node);
                        //cluster_sizes[max_block] += G.getNodeWeight(node);
                        //G.setPartitionIndex(node,max_block);
                //} endfor
        //}

        return 0;
}

void label_propagation_refinement::init_for_node_unit(graph_access& G, const uint64_t block_size,
                                                      const ParallelVector<Pair>& permutation,
                                                      Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                      std::unique_ptr<ConcurrentQueue>& queue) {
        Block block;
        block.reserve(block_size);
        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)].get().fetch_add(G.getNodeWeight(node),
                                                                         std::memory_order_relaxed);
                block.push_back(permutation[node].first);
                if (block.size() == block.capacity()) {
                        // block is full
                        queue->push(std::move(block));
                        block.clear();
                        block.reserve(block_size);
                }
        } endfor

        if (!block.empty()) {
                queue->push(std::move(block));
        }
}

template<typename T>
void label_propagation_refinement::seq_init_for_edge_unit(graph_access& G, const uint64_t block_size,
                                                      const T& permutation,
                                                      Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                      std::unique_ptr<ConcurrentQueue>& queue) {
        Block block;
        for (NodeID node = 0; node < G.number_of_nodes();) {
                uint32_t cur_block_size = 0;
                uint32_t size = 0;
                while (cur_block_size < block_size && node + size < G.number_of_nodes()) {
                        cur_block_size += permutation[node + size].second > 0 ? permutation[node + size].second : 1;
                        ++size;
                }
                block.reserve(size);
                NodeID end = node + size;
                for (; node < end; ++node) {
                        block.push_back(permutation[node].first);
                        cluster_sizes[G.getPartitionIndex(node)].get().fetch_add(G.getNodeWeight(node),
                                                                                 std::memory_order_relaxed);
                }

                queue->push(std::move(block));
                block.clear();
        }
        if (!block.empty()) {
                queue->push(std::move(block));
        }
}

template<typename T>
void label_propagation_refinement::par_init_for_edge_unit(graph_access& G, const uint64_t block_size,
                                                      const T& permutation,
                                                      Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                      std::unique_ptr<ConcurrentQueue>& queue) {

        CLOCK_START;
        std::atomic<NodeID> start_node(0);
        auto task = [&]() {
                const NodeID nodes_count = std::max<NodeID>(sqrt(G.number_of_nodes()), 4000);
                std::vector<NodeWeight> tmp_cluster_sizes(cluster_sizes.size());

                NodeID begin = start_node.fetch_add(nodes_count, std::memory_order_relaxed);

                while (begin < G.number_of_nodes()) {
                        NodeID end = std::min(begin + nodes_count, G.number_of_nodes());

                        for (auto node = begin; node != end; ++node) {
                                tmp_cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
                        }

                        begin = start_node.fetch_add(nodes_count, std::memory_order_relaxed);
                }

                return tmp_cluster_sizes;
        };

        std::vector<std::future<std::vector<NodeWeight>>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());
        for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                futures.emplace_back(parallel::g_thread_pool.Submit(i, task));
        }

        auto cur_cluster_sizes = task();
        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                auto tmp_cluster_sizes = future.get();
                for (size_t i = 0; i < cur_cluster_sizes.size(); ++i) {
                        cur_cluster_sizes[i] += tmp_cluster_sizes[i];
                }
        });

        for (size_t i = 0; i < cluster_sizes.size(); ++i) {
                cluster_sizes[i].get().store(cur_cluster_sizes[i], std::memory_order_relaxed);
        }
        CLOCK_END("Init cluster sizes");

        par_init_for_edge_unit(G, block_size, permutation, queue);
}

template<typename T>
void label_propagation_refinement::par_init_for_edge_unit(graph_access& G, const uint64_t block_size,
                                                          const T& permutation,
                                                          std::unique_ptr<ConcurrentQueue>& queue) {

        CLOCK_START;
        std::atomic<NodeID> offset(0);
        NodeID node_block_size = (NodeID) sqrt(G.number_of_nodes());
        node_block_size = std::max(node_block_size, 1000u);
        parallel::submit_for_all([&](uint32_t thread_id) {
                Block block;
                block.reserve(100);
                uint64_t cur_block_size = 0;
                while (true) {
                        NodeID begin = offset.fetch_add(node_block_size, std::memory_order_relaxed);
                        NodeID end = begin + node_block_size;
                        end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                        if (begin >= G.number_of_nodes()) {
                                break;
                        }

                        for (NodeID node = begin; node != end; ++node) {
                                block.push_back(permutation[node].first);
                                cur_block_size += permutation[node].second > 0 ? permutation[node].second : 1;
                                if (cur_block_size >= block_size) {
                                        // block is full
                                        queue->push(std::move(block));
                                        block.clear();
                                        block.reserve(100);
                                        cur_block_size = 0;
                                }
                        }
                }
                if (!block.empty()) {
                        queue->push(std::move(block));
                }
        });
        CLOCK_END("Init queue lp");
}


EdgeWeight label_propagation_refinement::parallel_label_propagation_with_queue(graph_access& G,
                                                                               PartitionConfig& config,
                                                                               Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                                               std::vector<std::vector<PartitionID>>& hash_maps,
                                                                               const parallel::ParallelVector<Pair>& permutation) {
        CLOCK_START;
        const NodeWeight block_upperbound = config.upper_bound_partition;
        auto queue = std::make_unique<ConcurrentQueue>();
        auto next_queue = std::make_unique<ConcurrentQueue>();

        std::vector<AtomicWrapper<bool>> queue_contains(G.number_of_nodes());
        std::vector<AtomicWrapper<bool>> next_queue_contains(G.number_of_nodes());

        uint32_t max_block_size = get_block_size(G, config);
        std::cout << "Uncoarsening: Block size\t" << max_block_size << std::endl;

//        TODO: check what is better with cache_aligned vector or not
//        std::vector<AtomicWrapper<NodeWeight>> cluster_sizes(config.k);

        bool use_edge_unit = config.block_size_unit == BlockSizeUnit::EDGES;
        if (use_edge_unit) {
                if (G.number_of_nodes() < 100) {
                        seq_init_for_edge_unit(G, max_block_size, permutation, cluster_sizes, queue);
                } else {
                        par_init_for_edge_unit(G, max_block_size, permutation, cluster_sizes, queue);
                }
        } else {
                init_for_node_unit(G, max_block_size, permutation, cluster_sizes, queue);
        }
        CLOCK_END("Uncoarsening: Parallel lp: Init queue lp");

        std::vector<std::future<NodeWeight>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());
        NodeWeight num_changed_label = 0;

        CLOCK_START_N;
        std::cout << "Uncoarsening: Num blocks\t" << queue->unsafe_size() << std::endl;
        for (int j = 0; j < config.label_iterations_refinement; j++) {
                if (queue->empty()) {
                        break;
                }
                auto process = [&](const size_t id) {
                        NodeWeight num_changed_label = 0;
                        Block cur_block;
                        Block new_block;
                        size_t new_block_size = 0;
                        new_block.reserve(100);

                        parallel::random rnd(config.seed + id);

                        std::vector<NodeID> neighbor_parts;

                        while (queue->try_pop(cur_block)) {
                                for (auto node : cur_block) {
                                        queue_contains[node].store(false, std::memory_order_relaxed);
                                        auto& hash_map = hash_maps[id];

                                        //now move the node to the cluster that is most common in the neighborhood
                                        neighbor_parts.clear();
                                        neighbor_parts.reserve(G.getNodeDegree(node));
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                PartitionID part = G.getPartitionIndex(target);
                                                if (hash_map[part] == 0) {
                                                        neighbor_parts.push_back(part);
                                                }
                                                hash_map[part] += G.getEdgeWeight(e);
                                        } endfor

                                        //second sweep for finding max and resetting array
                                        PartitionID my_block = G.getPartitionIndex(node);
                                        PartitionID max_block = my_block;
                                        NodeWeight max_cluster_size = 0;
                                        PartitionID max_value = 0;

                                        for (auto cur_part : neighbor_parts) {
                                                PartitionID cur_value = hash_map[cur_part];
                                                NodeWeight cur_cluster_size = cluster_sizes[cur_part].get().load(
                                                        std::memory_order_acquire);

                                                if ((cur_value > max_value || (cur_value == max_value && rnd.bit()))
                                                    &&
                                                    (cur_cluster_size + G.getNodeWeight(node) < block_upperbound
                                                     || (cur_part == my_block &&
                                                         cur_cluster_size <= block_upperbound))) {
                                                        max_value = cur_value;
                                                        max_block = cur_part;
                                                        max_cluster_size = cur_cluster_size;
                                                }

                                                hash_map[cur_part] = 0;
                                        }
                                        neighbor_parts.clear();


                                        bool changed_label = my_block != max_block;
                                        if (changed_label) {
                                                // try update size of the cluster
                                                bool perform_move = true;
                                                auto& atomic_val = cluster_sizes[max_block].get();
                                                while (!atomic_val.compare_exchange_weak(
                                                        max_cluster_size,
                                                        max_cluster_size + G.getNodeWeight(node),
                                                        std::memory_order_relaxed)) {
                                                        if (max_cluster_size + G.getNodeWeight(node) > block_upperbound) {
                                                                perform_move = false;
                                                                break;
                                                        }
                                                }

                                                if (perform_move) {
                                                        cluster_sizes[G.getPartitionIndex(node)].get().fetch_sub(
                                                                G.getNodeWeight(node), std::memory_order_relaxed);

                                                        G.setPartitionIndex(node, max_block);

                                                        ++num_changed_label;

                                                        forall_out_edges(G, e, node) {
                                                                NodeID target = G.getEdgeTarget(e);
                                                                if (!next_queue_contains[target].exchange(
                                                                        true, std::memory_order_relaxed)) {
                                                                        new_block.push_back(target);

                                                                        new_block_size += use_edge_unit ? G.getNodeDegree(target) : 1;
                                                                        if (new_block_size >= max_block_size) {
                                                                                next_queue->push(std::move(new_block));
                                                                                new_block.clear();
                                                                                new_block.reserve(100);
                                                                                new_block_size = 0;
                                                                        }
                                                                }
                                                        } endfor
                                                }
                                        }
                                }
                        }
                        if (!new_block.empty()) {
                                next_queue->push(std::move(new_block));
                        }
                        return num_changed_label;
                };


                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, process, i + 1));
                        //futures.push_back(parallel::g_thread_pool.Submit(process, i + 1));
                }

                //CLOCK_START;
                num_changed_label += process(0);
                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        num_changed_label += future.get();
                });
                //CLOCK_END(std::string("Iteration ") + std::to_string(j) + " time");
                std::swap(queue, next_queue);
                std::swap(queue_contains, next_queue_contains);
                futures.clear();
        }
        CLOCK_END("Uncoarsening: Parallel lp: iterations");

        return num_changed_label;
}

EdgeWeight label_propagation_refinement::parallel_label_propagation(graph_access& G,
                                                                    PartitionConfig& config,
                                                                    Cvector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                                    std::vector<std::vector<PartitionID>>& hash_maps,
                                                                    const parallel::ParallelVector<Pair>& permutation) {
        NodeWeight block_upperbound = config.upper_bound_partition;
        std::vector<std::future<NodeWeight>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());
        NodeWeight num_changed_label = 0;

        forall_nodes(G, node) {
                cluster_sizes[G.getPartitionIndex(node)].get().store(G.getNodeWeight(node), std::memory_order_relaxed);
        } endfor

        for (int j = 0; j < config.label_iterations_refinement; j++) {
                auto process = [&](const size_t id, NodeID begin, NodeID end) {
                        auto& hash_map = hash_maps[id];
                        NodeWeight num_changed_label = 0;
                        std::uniform_int_distribution<unsigned int> rnd(0,1);
                        std::mt19937 mt;
                        mt.seed(id);
                        for (NodeID node = begin; node != end; ++node) {
                                //now move the node to the cluster that is most common in the neighborhood
                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
                                } endfor

                                //second sweep for finding max and resetting array
                                PartitionID max_block = G.getPartitionIndex(node);
                                PartitionID my_block  = G.getPartitionIndex(node);
                                NodeWeight max_cluster_size = cluster_sizes[max_block].get().load(std::memory_order_acquire);

                                PartitionID max_value = 0;
                                forall_out_edges(G, e, node) {
                                        NodeID target             = G.getEdgeTarget(e);
                                        PartitionID cur_block     = G.getPartitionIndex(target);
                                        PartitionID cur_value     = hash_map[cur_block];
                                        NodeWeight cur_cluster_size = cluster_sizes[cur_block].get().load(std::memory_order_acquire);

                                        if((cur_value > max_value  || (cur_value == max_value
                                                                       && (bool) rnd(mt)))
                                           && (cur_cluster_size + G.getNodeWeight(node) < block_upperbound
                                               || (cur_block == my_block && cur_cluster_size <= block_upperbound)))
                                        {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                                max_cluster_size = cur_cluster_size;
                                        }

                                        hash_map[cur_block] = 0;
                                } endfor

                                bool changed_label = my_block != max_block;
                                if (changed_label) {
                                        // try update size of the cluster
                                        bool perform_move = true;
                                        auto& atomic_val = cluster_sizes[max_block].get();
                                        while (!atomic_val.compare_exchange_weak(
                                                max_cluster_size,
                                                max_cluster_size + G.getNodeWeight(node),
                                                std::memory_order_acq_rel)) {
                                                if (max_cluster_size + G.getNodeWeight(node) >
                                                    block_upperbound) {
                                                        perform_move = false;
                                                        break;
                                                }
                                        }

                                        if (perform_move) {
                                                cluster_sizes[G.getPartitionIndex(node)].get().fetch_sub(
                                                        G.getNodeWeight(node), std::memory_order_acq_rel);
                                                G.setPartitionIndex(node, max_block);
                                                ++num_changed_label;
                                        }
                                }
                        }
                        return num_changed_label;
                };

                size_t work_per_thread = G.number_of_nodes() / (parallel::g_thread_pool.NumThreads() + 1);
                NodeID first = 0;
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, process, i + 1, first,
                                                                         first + work_per_thread));
//                        futures.push_back(parallel::g_thread_pool.Submit(process, i + 1, first,
//                                                                         first + work_per_thread));

                        first += work_per_thread;
                }

                num_changed_label += process(0, first, G.number_of_nodes());
                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        num_changed_label += future.get();
                });
                futures.clear();
        }
        return num_changed_label;
}

EdgeWeight label_propagation_refinement::parallel_label_propagation_with_queue_with_many_clusters(graph_access& G,
                                                                               const PartitionConfig& config,
                                                                               const NodeWeight block_upperbound,
                                                                               std::vector<NodeWeight>& cluster_id,
                                                                               std::vector<AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                                               std::vector<std::vector<PartitionID>>& hash_maps,
                                                                               const parallel::ParallelVector<Triple>& permutation) {
        CLOCK_START;
        ALWAYS_ASSERT(config.block_size_unit == BlockSizeUnit::EDGES);
        auto queue = std::make_unique<ConcurrentQueue>();
        auto next_queue = std::make_unique<ConcurrentQueue>();

        std::vector<AtomicWrapper<bool>> queue_contains(G.number_of_nodes());
        std::vector<AtomicWrapper<bool>> next_queue_contains(G.number_of_nodes());

        uint64_t max_block_size = get_block_size(G, config);
        std::cout << "Block size\t" << max_block_size << std::endl;


        par_init_for_edge_unit(G, max_block_size, permutation, queue);
        CLOCK_END("Parallel lp: Init queue lp");

        NodeWeight num_changed_label = 0;

        ALWAYS_ASSERT(!config.graph_allready_partitioned);
        ALWAYS_ASSERT(!config.combine);

        using hash_function_type = parallel::MurmurHash<NodeID>;
        using hash_value_type = hash_function_type::hash_type;

        CLOCK_START_N;
        std::cout << "Num blocks\t" << queue->unsafe_size() << std::endl;
        for (int j = 0; j < config.label_iterations; j++) {
                if (queue->empty()) {
                        break;
                }

                auto process = [&](const size_t id) {
                        hash_function_type hash;
                        parallel::HashMap<NodeID, EdgeWeight, TabularHash<NodeID, 3, 2, 10, true>, true> hash_map(128);

                        NodeWeight num_changed_label = 0;
                        Block cur_block;
                        Block new_block;
                        size_t new_block_size = 0;
                        new_block.reserve(100);

                        parallel::random rnd(config.seed + id);
                        std::vector<NodeID> neighbor_parts;
                        while (queue->try_pop(cur_block)) {
                                for (auto node : cur_block) {
                                        hash.reset(config.seed + j + node);
                                        queue_contains[node].store(false, std::memory_order_relaxed);

                                        //now move the node to the cluster that is most common in the neighborhood
                                        neighbor_parts.clear();
                                        neighbor_parts.reserve(G.getNodeDegree(node));
                                        forall_out_edges(G, e, node){
                                                NodeID target = G.getEdgeTarget(e);
                                                NodeID cluster = cluster_id[target];
                                                auto& clst_size = hash_map[cluster];
                                                if (clst_size == 0) {
                                                        neighbor_parts.push_back(cluster);
                                                }
                                                clst_size += G.getEdgeWeight(e);
                                        } endfor

                                        //second sweep for finding max and resetting array
                                        const PartitionID my_block = cluster_id[node];
                                        PartitionID max_block = my_block;
                                        NodeWeight max_cluster_size = cluster_sizes[max_block].load(std::memory_order_relaxed);
                                        NodeWeight node_weight = G.getNodeWeight(node);
                                        EdgeWeight max_value = 0;
                                        hash_value_type max_block_hash = 0;
                                        for (auto cur_block : neighbor_parts) {
                                                EdgeWeight cur_value = hash_map[cur_block];
                                                NodeWeight cur_cluster_size = cluster_sizes[cur_block].load(std::memory_order_relaxed);
                                                hash_value_type cur_block_hash = 0;

                                                if ((cur_value > max_value || (cur_value == max_value && rnd.bit())) &&
                                                    (cur_cluster_size + node_weight < block_upperbound || cur_block == my_block)) {
                                                        max_value = cur_value;
                                                        max_block = cur_block;
                                                        max_cluster_size = cur_cluster_size;
                                                }
//                                                if ((cur_value > max_value || (cur_value == max_value && max_block_hash < (cur_block_hash = hash(cur_block)))) &&
//                                                    (cur_cluster_size + node_weight < block_upperbound || cur_block == my_block)) {
//                                                        if (cur_value > max_value) {
//                                                                cur_block_hash = hash(cur_block);
//                                                        }
//                                                        max_value = cur_value;
//                                                        max_block = cur_block;
//                                                        max_cluster_size = cur_cluster_size;
//                                                        max_block_hash = cur_block_hash;
//                                                }
                                        }
                                        hash_map.clear();

                                        bool changed_label = my_block != max_block;
                                        if (changed_label) {
                                                // try update size of the cluster
                                                bool perform_move = true;
                                                auto& atomic_val = cluster_sizes[max_block];
                                                while (!atomic_val.compare_exchange_weak(max_cluster_size,
                                                                                         max_cluster_size + node_weight,
                                                                                         std::memory_order_relaxed)) {
                                                        if (max_cluster_size + node_weight > block_upperbound) {
                                                                perform_move = false;
                                                                break;
                                                        }
                                                }

                                                if (perform_move) {
                                                        cluster_sizes[my_block].fetch_sub(node_weight,
                                                                                          std::memory_order_relaxed);

                                                        cluster_id[node] = max_block;

                                                        ++num_changed_label;

                                                        forall_out_edges(G, e, node) {
                                                                NodeID target = G.getEdgeTarget(e);
                                                                if (!next_queue_contains[target].exchange(
                                                                        true, std::memory_order_acq_rel)) {
                                                                        new_block.push_back(target);

                                                                        new_block_size += G.getNodeDegree(target);
                                                                        if (new_block_size >= max_block_size) {
                                                                                next_queue->push(std::move(new_block));
                                                                                new_block.clear();
                                                                                new_block.reserve(100);
                                                                                new_block_size = 0;
                                                                        }
                                                                }
                                                        } endfor
                                                }
                                        }
                                }
                        }
                        if (!new_block.empty()) {
                                next_queue->push(std::move(new_block));
                        }
                        return num_changed_label;
                };

                CLOCK_START;
                num_changed_label += parallel::submit_for_all(process, std::plus<NodeWeight>(), NodeWeight(0));
                CLOCK_END(std::string("Iteration ") + std::to_string(j) + " time");

                std::swap(queue, next_queue);
                std::swap(queue_contains, next_queue_contains);
        }
        CLOCK_END("Parallel lp: iterations");

        return num_changed_label;
}


EdgeWeight label_propagation_refinement::parallel_label_propagation_many_clusters(const PartitionConfig& config,
                                                                                  graph_access& G,
                                                                                  const NodeWeight block_upperbound,
                                                                                  std::vector<NodeWeight>& cluster_id,
                                                                                  NodeID& no_of_blocks) {
        CLOCK_START;
        std::vector<std::vector<PartitionID>> hash_maps(config.num_threads);
        std::vector<parallel::AtomicWrapper<NodeWeight>> cluster_sizes(G.number_of_nodes());

        forall_nodes(G, node) {
                cluster_id[node] = node;
                cluster_sizes[node].store(G.getNodeWeight(node), std::memory_order_relaxed);
        } endfor

        CLOCK_END("Init other vectors lp");

        CLOCK_START_N;
        parallel::ParallelVector<Triple> permutation(G.number_of_nodes());
        {
                CLOCK_START;
                std::atomic<NodeID> offset(0);
                NodeID block_size = (NodeID) sqrt(G.number_of_nodes());
                block_size = std::max(block_size, 1000u);

                parallel::submit_for_all([&](uint32_t thread_id) {
                        parallel::random rnd(config.seed + thread_id);
                        while (true) {
                                NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                                NodeID end = begin + block_size;
                                end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                                if (begin >= G.number_of_nodes()) {
                                        break;
                                }

                                for (NodeID node = begin; node != end; ++node) {
                                        permutation[node].first = node;
                                        permutation[node].second = G.getNodeDegree(node);
                                        permutation[node].rnd = rnd.random_number(0u, G.number_of_nodes() - 1);
                                }
                        }
                });

                CLOCK_END("Generate permutation");
        }

        parallel::g_thread_pool.Clear();
        parallel::Unpin();

        {
                CLOCK_START;
                parallel::sort(permutation.begin(), permutation.end(),
                               [](const Triple& lhs, const Triple& rhs) {
                                       return lhs.second < rhs.second
                                              || (lhs.second == rhs.second && lhs.rnd < rhs.rnd);
                               },
                               config.num_threads);
                CLOCK_END("Sort");
        }
        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(config.num_threads - 1);

        CLOCK_END("Parallel init of permutations lp");

        EdgeWeight res = 0;
        CLOCK_START_N;
        res = parallel_label_propagation_with_queue_with_many_clusters(G, config, block_upperbound, cluster_id,
                                                                       cluster_sizes, hash_maps, permutation);
        CLOCK_END("Main parallel (no queue) lp");

        std::cout << "Improved\t" << res << std::endl;

        CLOCK_START_N;
        if (config.num_threads > 1) {
                parallel_remap_cluster_ids_fast(config, G, cluster_id, no_of_blocks);
        } else {
                remap_cluster_ids_fast(config, G, cluster_id, no_of_blocks);
        }
        CLOCK_END("Remap cluster ids");
        return res;
}

void label_propagation_refinement::parallel_remap_cluster_ids_fast(const PartitionConfig& partition_config,
                                                                   graph_access& G,
                                                                   std::vector<NodeWeight>& cluster_id,
                                                                   NodeID& no_of_coarse_vertices) {
        if (cluster_id.empty()) {
                no_of_coarse_vertices = 0;
                return;
        }

        parallel::ParallelVector<AtomicWrapper<NodeID>> cluster_map(G.number_of_nodes());
        parallel::parallel_for_index(NodeID(0), G.number_of_nodes(), [&](NodeID node) {
                cluster_map[node] = 0;
        });

        parallel::parallel_for_index(NodeID(0), G.number_of_nodes(), [&](NodeID node) {
                PartitionID cur_cluster = cluster_id[node];
                cluster_map[cur_cluster] = 1;
        });

        parallel::partial_sum(cluster_map.begin(), cluster_map.end(), cluster_map.begin(),
                              partition_config.num_threads);

        parallel::parallel_for_index(NodeID(0), G.number_of_nodes(), [&](NodeID node) {
                cluster_id[node] = cluster_map[cluster_id[node]] - 1;
        });

        no_of_coarse_vertices = cluster_map.back();
}

void label_propagation_refinement::remap_cluster_ids_fast(const PartitionConfig& partition_config,
                                                          graph_access& G,
                                                          std::vector<NodeWeight>& cluster_id,
                                                          NodeID& no_of_coarse_vertices,
                                                          bool apply_to_graph) {
        if (cluster_id.empty()) {
                no_of_coarse_vertices = 0;
                return;
        }

        std::vector<NodeID> cluster_map(G.number_of_nodes());
        forall_nodes(G, node) {
                PartitionID cur_cluster = cluster_id[node];
                cluster_map[cur_cluster] = 1;
        } endfor

        for (size_t i = 1; i < cluster_map.size(); ++i)
                cluster_map[i] += cluster_map[i - 1];

        forall_nodes(G, node) {
                cluster_id[node] = cluster_map[cluster_id[node]] - 1;
        } endfor

        no_of_coarse_vertices = cluster_map.back();

        if (apply_to_graph) {
                forall_nodes(G, node) {
                                        G.setPartitionIndex(node, cluster_id[node]);
                                } endfor
                G.set_partition_count(no_of_coarse_vertices);
        }
}

EdgeWeight label_propagation_refinement::parallel_label_propagation(PartitionConfig& config, graph_access& G) {
        CLOCK_START;
        std::vector <std::vector<PartitionID>> hash_maps(config.num_threads, std::vector<PartitionID>(config.k));
        Cvector <AtomicWrapper<NodeWeight>> cluster_sizes(config.k);
        CLOCK_END("Uncoarsening: Init other vectors lp");

        CLOCK_START_N;
        parallel::ParallelVector<Pair> permutation(G.number_of_nodes());
        {
                CLOCK_START;
                std::atomic<NodeID> offset(0);
                NodeID block_size = (NodeID) sqrt(G.number_of_nodes());
                block_size = std::max(block_size, 1000u);

                parallel::submit_for_all([&](uint32_t thread_id) {
                        parallel::random rnd(config.seed + thread_id);
                        while (true) {
                                NodeID begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                                NodeID end = begin + block_size;
                                end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                                if (begin >= G.number_of_nodes()) {
                                        break;
                                }

                                for (NodeID node = begin; node != end; ++node) {
                                        permutation[node].first = node;
                                        permutation[node].second = G.getNodeDegree(node);
                                }
                        }
                });

                CLOCK_END("Generate permutation");
        }
        CLOCK_START_N;
        parallel::g_thread_pool.Clear();
        parallel::Unpin();
        {
                CLOCK_START;
                parallel::sort(permutation.begin(), permutation.end(),
                               [&](const Pair& lhs, const Pair& rhs) {
                                       return lhs.second < rhs.second
                                              || (lhs.second == rhs.second && lhs.first < rhs.first);
                               }, config.num_threads);
                CLOCK_END("Uncoarsening: Sort");
        }
        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(config.num_threads - 1);
        CLOCK_END("Uncoarsening: Sort with pool release");


        EdgeWeight res = 0;
        CLOCK_START_N;
        if (config.parallel_lp_type == ParallelLPType::NO_QUEUE) {
                res = parallel_label_propagation(G, config, cluster_sizes, hash_maps, permutation);
        } else if (config.parallel_lp_type == ParallelLPType::QUEUE) {
                res = parallel_label_propagation_with_queue(G, config, cluster_sizes, hash_maps, permutation);
        } else {
                res = 0;
        }
        CLOCK_END("Uncoarsening: Main parallel (no queue) lp");

        std::cout << "Uncoarsening: Improved\t" << res << std::endl;
        return res;
}

//void label_propagation_refinement::parallel_get_boundary_nodes(PartitionConfig& config, graph_access& G,
//                                                               std::vector<uint8_t>& boundary_nodes) {
//        std::vector<std::future<void>> futures;
//        futures.reserve(parallel::g_thread_pool.NumThreads());
//
//
//        uint32_t block_size = (uint32_t) sqrt(G.number_of_nodes());
//        block_size = std::max(block_size, 1000u);
//
//        std::atomic<uint32_t> offset(0);
//        auto process = [&](const uint32_t id) {
//                while (true) {
//                        size_t cur_index = offset.load(std::memory_order_relaxed);
//
//                        if (cur_index >= G.number_of_nodes()) {
//                                break;
//                        }
//
//                        uint32_t begin = offset.fetch_add(block_size, std::memory_order_relaxed);
//                        uint32_t end = begin + block_size;
//                        end = end <= G.number_of_nodes() ? end : G.number_of_nodes();
//
//                        if (begin >= G.number_of_nodes()) {
//                                break;
//                        }
//
//                        for (NodeID node = begin; node != end; ++node) {
//                                PartitionID cur_part = G.getPartitionIndex(node);
//                                forall_out_edges(G, e, node){
//                                        NodeID target = G.getEdgeTarget(e);
//                                        PartitionID part = G.getPartitionIndex(target);
//                                        if (cur_part != part) {
//                                                boundary_nodes[node] = (uint8_t) true;
//                                                break;
//                                        }
//                                } endfor
//                        }
//                }
//        };
//
//        for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
//                futures.push_back(parallel::g_thread_pool.Submit(i, process, i + 1));
//        }
//
//        process(0);
//        std::for_each(futures.begin(), futures.end(), [&](auto& future){
//                future.get();
//        });
//}
