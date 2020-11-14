#pragma once

#include <partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h>
#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"
#include "data_structure/parallel/random.h"

#include "matching.h"

namespace parallel {
class local_max_matching : public matching {
public:
        void match(const PartitionConfig& partition_config,
                   graph_access& G,
                   Matching& edge_matching,
                   CoarseMapping& mapping,
                   NodeID& no_of_coarse_vertices,
                   NodePermutationMap& permutation) override;

private:
        static constexpr NodeID m_none = std::numeric_limits<NodeID>::max();
        static constexpr uint32_t m_max_round = 20;

        using atomic_pair_type = std::pair<AtomicWrapper<NodeID>, AtomicWrapper<uint32_t>>;
        using block_type = std::vector<NodeID>;


        void parallel_match_with_queue(const PartitionConfig& partition_config,
                                       graph_access& G,
                                       Matching& edge_matching,
                                       CoarseMapping& mapping,
                                       NodeID& no_of_coarse_vertices);

        void parallel_match_with_queue_exp(const PartitionConfig& partition_config,
                                       graph_access& G,
                                       Matching& edge_matching,
                                       CoarseMapping& mapping,
                                       NodeID& no_of_coarse_vertices);

        void parallel_match_with_queue_exp_with_removal(const PartitionConfig& partition_config,
                                                        graph_access& G,
                                                        Matching& edge_matching,
                                                        CoarseMapping& mapping,
                                                        NodeID& no_of_coarse_vertices);

        void parallel_match(const PartitionConfig& partition_config,
                            graph_access& G,
                            Matching& edge_matching,
                            CoarseMapping& mapping,
                            NodeID& no_of_coarse_vertices);

        void sequential_match(const PartitionConfig& partition_config,
                              graph_access& G,
                              Matching& edge_matching,
                              CoarseMapping& mapping,
                              NodeID& no_of_coarse_vertices,
                              NodePermutationMap& permutation);

        void remap_matching(const PartitionConfig& partition_config, graph_access& G,
                            Matching& edge_matching, CoarseMapping& mapping, NodeID& no_of_coarse_vertices);

        NodeID find_max_neighbour_sequential(NodeID node, const uint32_t round, graph_access& G,
                                             const PartitionConfig& partition_config,
                                             std::vector<std::pair<NodeID, uint32_t>>& max_neighbours,
                                             Matching& vertex_mark, random& rnd) const;

        NodeID find_max_neighbour_parallel(NodeID node, const uint32_t round, graph_access& G,
                                           const PartitionConfig& partition_config,
                                           ParallelVector<atomic_pair_type>& max_neighbours,
                                           ParallelVector<AtomicWrapper<int>>& vertex_mark, random& rnd) const;

        NodeID find_max_neighbour_parallel(NodeID node, graph_access& G, const PartitionConfig& partition_config,
                                           ParallelVector<AtomicWrapper<int>>& vertex_mark, random& rnd) const;

        enum MatchingPhases {
                NOT_STARTED = 0,
                STARTED = 1,
                MATCHED = 2
        };
};
}