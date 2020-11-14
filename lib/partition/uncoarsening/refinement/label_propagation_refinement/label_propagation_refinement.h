/******************************************************************************
 * label_propagation_refinement.h  
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/


#ifndef LABEL_PROPAGATION_REFINEMENT_R4XW141Y
#define LABEL_PROPAGATION_REFINEMENT_R4XW141Y

#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"

#include "data_structure/parallel/pool_allocator.h"
#include "data_structure/parallel/thread_pool.h"

#include "definitions.h"
#include "../refinement.h"

#include <tbb/concurrent_queue.h>

#include <vector>

class label_propagation_refinement : public refinement {
public:
        label_propagation_refinement();
        virtual ~label_propagation_refinement();

        virtual EdgeWeight perform_refinement(PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary&);

        EdgeWeight perform_refinement(PartitionConfig & config, graph_access & G);

        EdgeWeight parallel_label_propagation_many_clusters(const PartitionConfig& config,
                                                            graph_access& G,
                                                            const NodeWeight block_upperbound,
                                                            std::vector<NodeWeight>& cluster_id,
                                                            NodeID& no_of_blocks);

private:
        using Block = std::vector<NodeID>;
        using ConcurrentQueue = tbb::concurrent_queue<Block>;
        using Pair = std::pair<NodeID, NodeID>;
        struct Triple {
                Triple(NodeID _first, NodeID _second, NodeID _rnd)
                        :       first(_first)
                        ,       second(_second)
                        ,       rnd(_rnd)
                {}

                NodeID first;
                NodeID second;
                NodeID rnd;
        };

        inline uint64_t get_block_size(graph_access& G, const PartitionConfig& config) const {
                if (config.block_size_unit == BlockSizeUnit::NODES) {
                        return get_block_size(G.number_of_nodes());
                }

                if (config.block_size_unit == BlockSizeUnit::EDGES) {
                        return get_block_size(G.number_of_edges());
                }

                return get_block_size(1);
        }

        inline uint64_t get_block_size(uint64_t total) const {
                uint64_t block_size = 1000;
                return std::max((uint64_t) sqrt(total), block_size);
        }

        std::chrono::system_clock::time_point begin, end;
        EdgeWeight sequential_label_propagation(PartitionConfig & config, graph_access & G);

        EdgeWeight parallel_label_propagation(PartitionConfig & config, graph_access & G);

        EdgeWeight parallel_label_propagation_with_queue(graph_access& G,
                                                         PartitionConfig& config,
                                                         parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                         std::vector<std::vector<PartitionID>>& hash_maps,
                                                         const parallel::ParallelVector<Pair>& permutation);

        EdgeWeight parallel_label_propagation_with_queue_with_many_clusters(graph_access& G,
                                                                            const PartitionConfig& config,
                                                                            const NodeWeight block_upperbound,
                                                                            std::vector<NodeWeight>& cluster_id,
                                                                            std::vector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                                                            std::vector<std::vector<PartitionID>>& hash_maps,
                                                                            const parallel::ParallelVector<Triple>& permutation);

        EdgeWeight parallel_label_propagation(graph_access& G,
                                              PartitionConfig& config,
                                              parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                              std::vector<std::vector<PartitionID>>& hash_maps,
                                              const parallel::ParallelVector<Pair>& permutation);

        template<typename T>
        void seq_init_for_edge_unit(graph_access& G, const uint64_t block_size,
                                    const T& permutation,
                                    parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                    std::unique_ptr<ConcurrentQueue>& queue);

        template<typename T>
        void par_init_for_edge_unit(graph_access& G, const uint64_t block_size,
                                    const T& permutation,
                                    parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                    std::unique_ptr<ConcurrentQueue>& queue);

        template<typename T>
        void par_init_for_edge_unit(graph_access& G, const uint64_t block_size,
                                    const T& permutation,
                                    std::unique_ptr<ConcurrentQueue>& queue);

        void init_for_node_unit(graph_access& G, const uint64_t block_size,
                                const parallel::ParallelVector<Pair>& permutation,
                                parallel::Cvector<parallel::AtomicWrapper<NodeWeight>>& cluster_sizes,
                                std::unique_ptr<ConcurrentQueue>& queue);

//        void parallel_get_boundary_nodes(PartitionConfig& config, graph_access& G,
//                                         std::vector<uint8_t>& bounday_nodes);

        void remap_cluster_ids_fast(const PartitionConfig& partition_config, graph_access& G,
                                    std::vector<NodeWeight>& cluster_id, NodeID& no_of_coarse_vertices,
                                    bool apply_to_graph = false);

        void parallel_remap_cluster_ids_fast(const PartitionConfig& partition_config, graph_access& G,
                                             std::vector<NodeWeight>& cluster_id, NodeID& no_of_coarse_vertices);
};


#endif /* end of include guard: LABEL_PROPAGATION_REFINEMENT_R4XW141Y */
