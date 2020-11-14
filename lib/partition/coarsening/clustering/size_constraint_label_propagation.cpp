/******************************************************************************
 * size_constraint_label_propagation.h     
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


#include <unordered_map>

#include <sstream>
#include "../edge_rating/edge_ratings.h"
#include "../matching/gpa/gpa_matching.h"
#include "data_structure/union_find.h"
#include "data_structure/parallel/time.h"
#include "node_ordering.h"
#include "partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.h"
#include "partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "tools/quality_metrics.h"
#include "tools/random_functions.h"
#include "io/graph_io.h"

#include "size_constraint_label_propagation.h"

#include <ips4o/ips4o.hpp>

//#include <parallel/algorithm>

size_constraint_label_propagation::size_constraint_label_propagation() {
                
}

size_constraint_label_propagation::~size_constraint_label_propagation() {
                
}

void size_constraint_label_propagation::match(const PartitionConfig & partition_config, 
                                              graph_access & G, 
                                              Matching & _matching, 
                                              CoarseMapping & coarse_mapping, 
                                              NodeID & no_of_coarse_vertices,
                                              NodePermutationMap & permutation) {
        coarse_mapping.resize(G.number_of_nodes());
        no_of_coarse_vertices = 0;

        if ( partition_config.ensemble_clusterings ) {
                ensemble_clusterings(partition_config, G, _matching, coarse_mapping, no_of_coarse_vertices, permutation);
        } else {
                match_internal(partition_config, G, _matching, coarse_mapping, no_of_coarse_vertices, permutation);
        }
}

void size_constraint_label_propagation::match_internal(const PartitionConfig & partition_config, 
                                              graph_access & G,
                                              Matching & _matching,
                                              CoarseMapping & coarse_mapping,
                                              NodeID & no_of_coarse_vertices,
                                              NodePermutationMap&) {

        std::vector<NodeWeight> cluster_id(G.number_of_nodes());
        NodeWeight block_upperbound = ceil(partition_config.upper_bound_partition/(double)partition_config.cluster_coarsening_factor);
        std::cout << "BLOCK UPPER BOUND = " << block_upperbound << std::endl;
        if (!partition_config.parallel_coarsening_lp) {
                label_propagation(partition_config, G, block_upperbound, cluster_id, no_of_coarse_vertices);
        } else
        {
//                parallel_label_propagation(partition_config, G, block_upperbound, cluster_id, no_of_coarse_vertices);
                label_propagation_refinement().parallel_label_propagation_many_clusters(partition_config, G,
                                                                                        block_upperbound, cluster_id, no_of_coarse_vertices);
        }
        create_coarsemapping( partition_config, G, cluster_id, coarse_mapping);
}

void size_constraint_label_propagation::ensemble_two_clusterings( graph_access & G, 
                                                                  std::vector<NodeID> & lhs, 
                                                                  std::vector<NodeID> & rhs, 
                                                                  std::vector< NodeID > & output,
                                                                  NodeID & no_of_coarse_vertices) {


        hash_ensemble new_mapping; 
        no_of_coarse_vertices = 0;
        for( NodeID node = 0; node < lhs.size(); node++) {
                ensemble_pair cur_pair;
                cur_pair.lhs = lhs[node]; 
                cur_pair.rhs = rhs[node]; 
                cur_pair.n   = G.number_of_nodes(); 

                if(new_mapping.find(cur_pair) == new_mapping.end() ) {
                        new_mapping[cur_pair].mapping = no_of_coarse_vertices;
                        no_of_coarse_vertices++;
                }

                output[node] = new_mapping[cur_pair].mapping;
        }

        no_of_coarse_vertices = new_mapping.size();
}


void size_constraint_label_propagation::ensemble_clusterings(const PartitionConfig & partition_config, 
                                                             graph_access & G, 
                                                             Matching & _matching, 
                                                             CoarseMapping & coarse_mapping, 
                                                             NodeID & no_of_coarse_vertices,
                                                             NodePermutationMap &) {
        int runs = partition_config.number_of_clusterings;
        std::vector< NodeID >  cur_cluster(G.number_of_nodes(), 0);
        std::vector< NodeID >  ensemble_cluster(G.number_of_nodes(),0);

        int new_cf = partition_config.cluster_coarsening_factor;
        for( int i = 0; i < runs; i++) {
                PartitionConfig config = partition_config;
                config.cluster_coarsening_factor = new_cf;

                NodeID cur_no_blocks = 0;
                label_propagation(config, G, cur_cluster, cur_no_blocks); 

                if( i != 0 ) {
                        ensemble_two_clusterings(G, cur_cluster, ensemble_cluster, ensemble_cluster, no_of_coarse_vertices);
                } else {
                        forall_nodes(G, node) {
                                ensemble_cluster[node] = cur_cluster[node];
                        } endfor
                        
                        no_of_coarse_vertices = cur_no_blocks;
                }
                new_cf = random_functions::nextInt(10, 30);
        }

        create_coarsemapping( partition_config, G, ensemble_cluster, coarse_mapping);


}

void size_constraint_label_propagation::label_propagation(const PartitionConfig & partition_config, 
                                                         graph_access & G, 
                                                         std::vector<NodeWeight> & cluster_id, 
                                                         NodeID & no_of_blocks ) {
        NodeWeight block_upperbound = ceil(partition_config.upper_bound_partition/(double)partition_config.cluster_coarsening_factor);

        label_propagation( partition_config, G, block_upperbound, cluster_id, no_of_blocks);
}

void size_constraint_label_propagation::label_propagation(const PartitionConfig & partition_config, 
                                                         graph_access & G, 
                                                         const NodeWeight & block_upperbound,
                                                         std::vector<NodeWeight> & cluster_id,  
                                                         NodeID & no_of_blocks) {
        // in this case the _matching paramter is not used 
        // coarse_mappng stores cluster id and the mapping (it is identical)
        std::vector<PartitionID> hash_map(G.number_of_nodes(),0);
        std::vector<NodeID> permutation(G.number_of_nodes());
        std::vector<NodeWeight> cluster_sizes(G.number_of_nodes());
        cluster_id.resize(G.number_of_nodes());

        forall_nodes(G, node) {
                cluster_sizes[node] = G.getNodeWeight(node);
                cluster_id[node]    = node;
        } endfor
        
        node_ordering n_ordering;
        n_ordering.order_nodes(partition_config, G, permutation);

        for( int j = 0; j < partition_config.label_iterations; j++) {
                unsigned int change_counter = 0;
                forall_nodes(G, i) {
                        NodeID node = permutation[i];
                        //now move the node to the cluster that is most common in the neighborhood

                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                hash_map[cluster_id[target]]+=G.getEdgeWeight(e);
                        } endfor

                        //second sweep for finding max and resetting array
                        PartitionID max_block = cluster_id[node];
                        PartitionID my_block  = cluster_id[node];

                        PartitionID max_value = 0;
                        forall_out_edges(G, e, node) {
                                NodeID target             = G.getEdgeTarget(e);
                                PartitionID cur_block     = cluster_id[target];
                                PartitionID cur_value     = hash_map[cur_block];
                                if((cur_value > max_value  || (cur_value == max_value && random_functions::nextBool())) 
                                && (cluster_sizes[cur_block] + G.getNodeWeight(node) < block_upperbound || cur_block == my_block) 
                                && (!partition_config.graph_allready_partitioned || G.getPartitionIndex(node) == G.getPartitionIndex(target))
                                && (!partition_config.combine || G.getSecondPartitionIndex(node) == G.getSecondPartitionIndex(target)))
                                {
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }

                                hash_map[cur_block] = 0;
                        } endfor

                        cluster_sizes[cluster_id[node]]  -= G.getNodeWeight(node);
                        cluster_sizes[max_block]         += G.getNodeWeight(node);
                        change_counter                   += (cluster_id[node] != max_block);
                        cluster_id[node]                  = max_block;
                } endfor
        }

        remap_cluster_ids_fast(partition_config, G, cluster_id, no_of_blocks);
}

uint32_t size_constraint_label_propagation::parallel_label_propagation(const PartitionConfig& config,
                                                                       graph_access& G,
                                                                       const NodeWeight block_upperbound,
                                                                       std::vector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                                       std::vector<NodeID>& cluster_id,
                                                                       std::vector<pair_type>& permutation,
                                                                       NodeID& no_of_blocks,
                                                                       std::vector<parallel::AtomicWrapper<char>> active,
                                                                       std::vector<parallel::AtomicWrapper<char>> new_active
) {
        std::vector<std::future<std::pair<uint32_t, uint32_t>>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());
        uint32_t num_changed_label_all = 0;

        std::vector<std::vector<PartitionID>> hash_maps(config.num_threads);

        uint32_t block_size = (uint32_t) sqrt(G.number_of_nodes());
        block_size = std::max(block_size, 1000u);

        for (int j = 0; j < config.label_iterations; j++) {
                std::atomic<uint32_t> offset(0);
                auto process = [&](const uint32_t id) {
                        uint32_t num_active = 0;
                        auto& hash_map = hash_maps[id];

                        if (hash_map.empty()) {
                                hash_map.resize(G.number_of_nodes());
                        }

                        uint32_t num_changed_label = 0;
                        parallel::random rnd(config.seed + id);
                        std::vector<NodeID> neighbor_parts;

                        while (true) {
                                size_t cur_index = offset.load(std::memory_order_relaxed);

                                if (cur_index >= G.number_of_nodes()) {
                                        break;
                                }

                                uint32_t begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                                uint32_t end = begin + block_size;
                                end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                                if (begin >= G.number_of_nodes()) {
                                        break;
                                }

                                for (NodeID index = begin; index != end; ++index) {
                                        NodeID node = permutation[index].first;

                                        if (!active[node].load(std::memory_order_relaxed)) {
                                                continue;
                                        }
                                        active[node].store(false, std::memory_order_relaxed);

                                        const PartitionID my_block = cluster_id[node];
                                        //now move the node to the cluster that is most common in the neighborhood
                                        neighbor_parts.clear();
                                        neighbor_parts.reserve(G.getNodeDegree(node));
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                NodeID cluster = cluster_id[target];
                                                if (hash_map[cluster] == 0) {
                                                        neighbor_parts.push_back(cluster);
                                                }
                                                hash_map[cluster] += G.getEdgeWeight(e);
                                        } endfor

                                        //second sweep for finding max and resetting array
                                        PartitionID max_block = my_block;
                                        NodeWeight max_cluster_size = cluster_sizes[max_block].load(
                                                std::memory_order_relaxed);

                                        PartitionID max_value = 0;
                                        NodeWeight node_weight = G.getNodeWeight(node);
                                        for (auto cur_block : neighbor_parts) {
                                                PartitionID cur_value = hash_map[cur_block];

                                                NodeWeight cur_cluster_size = cluster_sizes[cur_block].load(
                                                        std::memory_order_relaxed);

                                                if ((cur_value > max_value || (cur_value == max_value && rnd.bit())) &&
                                                        (cur_cluster_size + node_weight < block_upperbound || cur_block == my_block)) {
                                                        ALWAYS_ASSERT(!config.graph_allready_partitioned);
                                                        ALWAYS_ASSERT(!config.combine);
                                                        max_value = cur_value;
                                                        max_block = cur_block;
                                                        max_cluster_size = cur_cluster_size;
                                                }

                                                hash_map[cur_block] = 0;
                                        }

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
                                                                new_active[target].store(true, std::memory_order_relaxed);
//                                                                char res = false;
//                                                                if (new_active[target].compare_exchange_strong(res, true, std::memory_order_relaxed)) {
//                                                                        ++num_active;
//                                                                }
                                                        } endfor
                                                }
                                        }
                                }
                        }
                        return std::make_pair(num_changed_label, num_active);
                };

                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, process, i + 1));
//                        futures.push_back(parallel::g_thread_pool.Submit(process, i + 1));
                }

                uint32_t num_active = 0;
                uint32_t na = 0;
                uint32_t changed = 0;

                CLOCK_START;
                std::tie(changed, na) = process(0);
                num_changed_label_all += changed;
                num_active += na;

                std::for_each(futures.begin(), futures.end(), [&](auto& future){
                        std::tie(changed, na) = future.get();
                        num_changed_label_all += changed;
                        num_active += na;
                });
                CLOCK_END(std::string("Iteration ") + std::to_string(j) + " time");
                active.swap(new_active);
                block_size = std::max((uint32_t) sqrt(num_active), 1000u);
                futures.clear();
        }
        return num_changed_label_all;
}

//uint32_t size_constraint_label_propagation::parallel_label_propagation_exp(const PartitionConfig& config,
//                                                                       graph_access& G,
//                                                                       const NodeWeight block_upperbound,
//                                                                       std::vector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
//                                                                       std::vector<NodeID>& cluster_id,
//                                                                       std::vector<pair_type>& permutation,
//                                                                       NodeID& no_of_blocks,
//                                                                       std::vector<parallel::AtomicWrapper<char>> active,
//                                                                       std::vector<parallel::AtomicWrapper<char>> new_active
//) {
//        std::vector<std::unique_ptr<PartitionID[]>> hash_maps(config.num_threads);
//        std::vector<std::vector<NodeID>> neighbor_parts_array(config.num_threads);
//        parallel::Cvector<parallel::random> rnds;
//        rnds.reserve(config.num_threads);
//
//        for (size_t tid = 0; tid < config.num_threads; ++tid) {
//                rnds.emplace_back(config.seed + tid);
//        }
//
//        for (int j = 0; j < config.label_iterations; j++) {
//                auto task = [&](const auto& node_degree_pair) {
//                        int id = omp_get_thread_num();
//                        auto& hash_map = hash_maps[id];
//                        auto& neighbor_parts = neighbor_parts_array[id];
//                        auto& rnd = rnds[id].get();
//
//                        if (hash_map.get() == nullptr) {
//                                hash_map = std::make_unique<PartitionID[]>(G.number_of_nodes());
//                        }
//
//                        NodeID node = node_degree_pair.first;
//
//                        if (!active[node].load(std::memory_order_relaxed)) {
//                                return;
//                        }
//                        active[node].store(false, std::memory_order_relaxed);
//
//                        const PartitionID my_block = cluster_id[node];
//                        //now move the node to the cluster that is most common in the neighborhood
//                        neighbor_parts.clear();
//                        neighbor_parts.reserve(std::max<uint32_t>(G.getNodeDegree(node), parallel::g_cache_line_size / sizeof(NodeID)));
//                        forall_out_edges(G, e, node) {
//                                NodeID target = G.getEdgeTarget(e);
//                                NodeID cluster = cluster_id[target];
//                                if (hash_map[cluster] == 0) {
//                                        neighbor_parts.push_back(cluster);
//                                }
//                                hash_map[cluster] += G.getEdgeWeight(e);
//                        } endfor
//
//                        //second sweep for finding max and resetting array
//                        PartitionID max_block = my_block;
//                        NodeWeight max_cluster_size = cluster_sizes[max_block].load(
//                                std::memory_order_relaxed);
//
//                        PartitionID max_value = 0;
//                        NodeWeight node_weight = G.getNodeWeight(node);
//                        for (auto cur_block : neighbor_parts) {
//                                PartitionID cur_value = hash_map[cur_block];
//
//                                NodeWeight cur_cluster_size = cluster_sizes[cur_block].load(
//                                        std::memory_order_relaxed);
//
//                                if ((cur_value > max_value || (cur_value == max_value && rnd.bit())) &&
//                                    (cur_cluster_size + node_weight < block_upperbound || cur_block == my_block)) {
//                                        ALWAYS_ASSERT(!config.graph_allready_partitioned);
//                                        ALWAYS_ASSERT(!config.combine);
//                                        max_value = cur_value;
//                                        max_block = cur_block;
//                                        max_cluster_size = cur_cluster_size;
//                                }
//
//                                hash_map[cur_block] = 0;
//                        }
//
//                        bool changed_label = my_block != max_block;
//                        if (changed_label) {
//                                // try update size of the cluster
//                                bool perform_move = true;
//                                auto& atomic_val = cluster_sizes[max_block];
//                                while (!atomic_val.compare_exchange_weak(max_cluster_size,
//                                                                         max_cluster_size + node_weight,
//                                                                         std::memory_order_relaxed)) {
//                                        if (max_cluster_size + node_weight > block_upperbound) {
//                                                perform_move = false;
//                                                break;
//                                        }
//                                }
//
//                                if (perform_move) {
//                                        cluster_sizes[my_block].fetch_sub(node_weight,
//                                                                          std::memory_order_relaxed);
//
//                                        cluster_id[node] = max_block;
//
//                                        forall_out_edges(G, e, node) {
//                                                NodeID target = G.getEdgeTarget(e);
//                                                new_active[target].store(true, std::memory_order_relaxed);
//                                        } endfor
//                                }
//                        }
//                };
//
//                CLOCK_START;
//                __gnu_parallel::for_each(permutation.begin(), permutation.end(), task, __gnu_parallel::_Parallelism::parallel_balanced);
//                CLOCK_END(std::string("Iteration ") + std::to_string(j) + " time");
//                active.swap(new_active);
//        }
//        return 0;
//}

void size_constraint_label_propagation::parallel_label_propagation(const PartitionConfig& config,
                                                                   graph_access& G,
                                                                   const NodeWeight block_upperbound,
                                                                   std::vector<NodeWeight>& cluster_id,
                                                                   NodeID& no_of_blocks) {
        CLOCK_START;
        std::vector<parallel::AtomicWrapper<NodeWeight>> cluster_sizes(G.number_of_nodes());
        std::vector<parallel::AtomicWrapper<char>> active(G.number_of_nodes(), true);
        std::vector<parallel::AtomicWrapper<char>> new_active(G.number_of_nodes(), false);

        forall_nodes(G, node) {
                cluster_id[node] = node;
                cluster_sizes[node].store(G.getNodeWeight(node), std::memory_order_relaxed);
        } endfor

        CLOCK_END("Init other vectors lp");

        CLOCK_START_N;
        std::vector<pair_type> permutation;
        permutation.reserve(G.number_of_nodes());

        forall_nodes(G, node) {
                permutation.emplace_back(node, G.getNodeDegree(node));
        } endfor

        parallel::g_thread_pool.Clear();
        parallel::Unpin();
        {
                CLOCK_START;
                parallel::sort(permutation.begin(), permutation.end(), [&](const pair_type& lhs, const pair_type& rhs) {
                        return lhs.second < rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
                }, config.num_threads);
                CLOCK_END("Sort");
        }
        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(config.num_threads - 1);

        {
                CLOCK_START;
                parallel::random rnd(config.seed);
                size_t i = 0;
                while (i < permutation.size()) {
                        size_t end = i;
                        while (end + 1 < permutation.size() && permutation[end].second == permutation[end + 1].second) {
                                ++end;
                        }

                        rnd.shuffle(permutation.begin() + i, permutation.begin() + end + 1);
                        i = end + 1;
                }
                CLOCK_END("Shuffle");
        }

        CLOCK_END("Parallel init of permutations lp");

        CLOCK_START_N;

        uint32_t num_changed_label = parallel_label_propagation(config, G, block_upperbound, cluster_sizes, cluster_id,
                                                                permutation, no_of_blocks, active, new_active);

        CLOCK_END("Main parallel (no queue) lp");
        std::cout << "Improved\t" << num_changed_label << std::endl;

        CLOCK_START_N;
        if (config.num_threads > 1) {
                parallel_remap_cluster_ids_fast(config, G, cluster_id, no_of_blocks);
        } else {
                remap_cluster_ids_fast(config, G, cluster_id, no_of_blocks);
        }
        CLOCK_END("Remap cluster ids");
}

void size_constraint_label_propagation::create_coarsemapping(const PartitionConfig & partition_config, 
                                                             graph_access & G,
                                                             std::vector<NodeWeight>& cluster_id,
                                                             CoarseMapping & coarse_mapping) {
        forall_nodes(G, node) {
                coarse_mapping[node] = cluster_id[node];
        } endfor
}

void size_constraint_label_propagation::remap_cluster_ids(const PartitionConfig & partition_config, 
                                                          graph_access & G,
                                                          std::vector<NodeWeight> & cluster_id,
                                                          NodeID & no_of_coarse_vertices, bool apply_to_graph) {

        PartitionID cur_no_clusters = 0;
        std::unordered_map<PartitionID, PartitionID> remap;
        forall_nodes(G, node) {
                PartitionID cur_cluster = cluster_id[node];
                //check wether we already had that
                if( remap.find( cur_cluster ) == remap.end() ) {
                        remap[cur_cluster] = cur_no_clusters++;
                }

                cluster_id[node] = remap[cur_cluster];
        } endfor

        if( apply_to_graph ) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, cluster_id[node]);
                } endfor
                G.set_partition_count(cur_no_clusters);
        }

        no_of_coarse_vertices = cur_no_clusters;
}

void size_constraint_label_propagation::remap_cluster_ids_fast(const PartitionConfig& partition_config,
                                                               graph_access& G,
                                                               std::vector<NodeWeight>& cluster_id,
                                                               NodeID& no_of_coarse_vertices,
                                                               bool apply_to_graph) {
        if (cluster_id.empty()) {
                no_of_coarse_vertices = 0;
                return;
        }

        std::vector<NodeWeight> cluster_map(G.number_of_nodes());
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

void size_constraint_label_propagation::parallel_remap_cluster_ids_fast(const PartitionConfig& partition_config,
                                                                        graph_access& G,
                                                                        std::vector<NodeWeight>& cluster_id,
                                                                        NodeID& no_of_coarse_vertices) {
        if (cluster_id.empty()) {
                no_of_coarse_vertices = 0;
                return;
        }

        parallel::ParallelVector<NodeID> cluster_map(G.number_of_nodes());
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
