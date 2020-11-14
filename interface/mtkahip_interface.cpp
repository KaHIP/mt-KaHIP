/******************************************************************************
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

#include <iostream>
#include "mtkahip_interface.h"
#include "../lib/data_structure/graph_access.h"
#include "../lib/io/graph_io.h"
#include "../lib/tools/timer.h"
#include "../lib/tools/quality_metrics.h"
#include "../lib/tools/macros_assertions.h"
#include "../lib/tools/random_functions.h"
#include "../lib/parallel_mh/parallel_mh_async.h"
#include "../lib/partition/uncoarsening/separator/area_bfs.h"
#include "../lib/partition/partition_config.h"
#include "../lib/partition/graph_partitioner.h"
#include "../lib/partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "../app/configuration.h"
#include "../app/balance_configuration.h"
#include "../data_structure/parallel/thread_pool.h"

#ifdef __gnu_linux__
#include <numa.h>
#endif

using namespace std;

void internal_build_graph( PartitionConfig & partition_config, 
                           int* n, 
                           int* vwgt, 
                           int* xadj, 
                           int* adjcwgt, 
                           int* adjncy,
                           graph_access & G) {
        G.build_from_metis(*n, xadj, adjncy); 
        G.set_partition_count(partition_config.k); 
 
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
       
        if(vwgt != NULL) {
                forall_nodes(G, node) {
                        G.setNodeWeight(node, vwgt[node]);
                } endfor
        }

        if(adjcwgt != NULL) {
                forall_edges(G, e) {
                        G.setEdgeWeight(e, adjcwgt[e]);
                } endfor 
        }

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);
}

void internal_mtkahip_call(PartitionConfig & partition_config, 
                          bool suppress_output, 
                          int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int* nparts, 
                          double* imbalance, 
                          int* edgecut, 
                          int* part) {

        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        //if(suppress_output) {
        cout.rdbuf(ofs.rdbuf()); 
        //}

        partition_config.imbalance = 100*(*imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

        graph_partitioner partitioner;
        partitioner.perform_partitioning(partition_config, G);

        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);

        ofs.close();
        cout.rdbuf(backup);
}

void mtkahip(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool suppress_output, 
                   int seed,
                   int mode,
                   uint32_t num_threads,
                   int* edgecut, 
                   int* part) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;
        partition_config.num_threads = num_threads;

#ifdef __gnu_linux__
        bitmask* old_numa_mask_ptr = nullptr;
#endif

#ifdef __gnu_linux__
        if (numa_available() < 0) {
                printf("No NUMA support available on this system.\n");
                exit(1);
        }
        numa_set_interleave_mask(numa_all_nodes_ptr);
        old_numa_mask_ptr = numa_get_interleave_mask();
#endif
        if (mode == FASTSOCIALMULTITRY_PARALLEL) {
                cfg.fastsocialmultitry_parallel(partition_config);
        } else if (mode == FASTSOCIAL_PARALLEL) {
                cfg.fastsocial_parallel(partition_config);
        }

        parallel::PinToCore(partition_config.main_core);
        parallel::g_thread_pool.Resize(partition_config.num_threads - 1);

        partition_config.seed = seed;
        internal_mtkahip_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, edgecut, part);

        if (mode == FASTSOCIALMULTITRY_PARALLEL || mode == FASTSOCIAL_PARALLEL) {
                parallel::Unpin();
                parallel::g_thread_pool.Clear();
#ifdef __gnu_linux__
                numa_set_interleave_mask(old_numa_mask_ptr);
#endif
        }
}



