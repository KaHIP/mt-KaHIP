#pragma once

#include "data_structure/graph_access.h"
#include "partition_config.h"


namespace parallel {

class initial_partitioning {
public:
        initial_partitioning();

        virtual ~initial_partitioning();

        void perform_initial_partitioning(const PartitionConfig& config, graph_access& G);
};

}