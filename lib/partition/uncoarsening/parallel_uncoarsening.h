#include "data_structure/graph_hierarchy.h"
#include "partition/partition_config.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/fast_boundary.h"

namespace parallel {

class uncoarsening {
public:
        int perform_uncoarsening_cut(const PartitionConfig& config, graph_hierarchy& hierarchy);
private:
        void perform_label_propagation(PartitionConfig& config, graph_access& G);

        EdgeWeight perform_multitry_kway(PartitionConfig& config, graph_access& G, boundary_type& boundary);
};

}