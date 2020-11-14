#pragma once

#include "data_structure/graph_access.h"
#include "partition/coarsening/matching/matching.h"

class hash_common_neighborhood : public matching {
public:
        void match(const PartitionConfig& config,
                   graph_access& G,
                   Matching&,
                   CoarseMapping& coarse_mapping,
                   NodeID& no_of_coarse_vertices,
                   NodePermutationMap&) override;

        virtual ~hash_common_neighborhood() {}

private:
        void find_vertices_with_common_neighbors(const PartitionConfig& config,
                                                 graph_access& G,
                                                 CoarseMapping& coarse_mapping,
                                                 NodeID& no_of_coarse_vertices);
};
