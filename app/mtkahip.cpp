/******************************************************************************
 * kaffpa.cpp 
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

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/parallel/graph_utils.h"
#include "data_structure/parallel/thread_pool.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_stop_rule.h"
#include "uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/multitry_kway_fm.h"
#include "uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_commons.h"

#include "data_structure/parallel/adaptive_hash_table.h"

#ifdef __gnu_linux__
#include <numa.h>
#endif

#include <execinfo.h>
#include <omp.h>
#include <signal.h>

void handler(int sig) {
        void *array[10];
        size_t size;

        // get void*'s for all entries on the stack
        size = backtrace(array, 10);

        // print out all the frames to stderr
        fprintf(stderr, "Error: signal %d:\n", sig);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
        exit(1);
}

int main(int argn, char **argv) {
        omp_set_dynamic(false);
        omp_set_num_threads(0);

        signal(SIGSEGV, handler);
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
        #pragma message("COMPARE WITH SEQUENTIAL MODE IS ON")
        std::cout << "COMPARE WITH SEQUENTIAL MODE IS ON!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << "THIS SLOWS DOWN THE APP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::end;
#endif

#ifdef __gnu_linux__
        if (numa_available() < 0) {
		printf("No NUMA support available on this system.\n");
		exit(1);
	}
        numa_set_interleave_mask(numa_all_nodes_ptr);
#endif

        //std::cout << "Git revision\t" << GIT_DESC << std::endl;
        PartitionConfig partition_config;
        std::string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;
       
        int ret_code = parse_parameters(argn, argv, 
                                        partition_config, 
                                        graph_filename, 
                                        is_graph_weighted, 
                                        suppress_output, recursive);

        if(ret_code) {
                return 0;
        }

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);
        graph_access G;     
        ALWAYS_ASSERT(partition_config.main_core == 0);

        timer t;
        if (!partition_config.shuffle_graph && !partition_config.sort_edges) {
                graph_io::readGraphWeighted(G, graph_filename);
                //double avg;
                //double med;
                //std::tie(avg, med) = average_chain_length(G);
                //std::cout << "avg median chain len = " << avg << std::endl;
                //std::cout << "median median chain len = " << med << std::endl;
        } else {
                ALWAYS_ASSERT(!partition_config.shuffle_graph || !partition_config.sort_edges);
                if (partition_config.shuffle_graph) {
                        graph_access tmp_G;
                        graph_io::readGraphWeighted(tmp_G, graph_filename);
                        shuffle_graph(tmp_G, G);
                }
                if (partition_config.sort_edges) {
                        graph_access tmp_G;
                        graph_io::readGraphWeighted(tmp_G, graph_filename);
                        sort_edges(tmp_G, G);
                }
        }
        std::cout << "io time: " << t.elapsed()  << std::endl;
        G.set_partition_count(partition_config.k);
 
        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        quality_metrics qm;
        Gain input_partition_cut = 0;
        std::vector<PartitionID> input_partition;
        if (partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned = true;
                partition_config.only_first_level = true;
                partition_config.mh_no_mh = false;
                partition_config.no_change_convergence = false;
                partition_config.corner_refinement_enabled = false;
                partition_config.kaffpa_perfectly_balanced_refinement = false;
                partition_config.quotient_graph_two_way_refinement = false;

                input_partition.resize(G.number_of_nodes());

                forall_nodes(G, node) {
                        input_partition[node] = G.getPartitionIndex(node);
                } endfor
                input_partition_cut = qm.edge_cut(G);
        }
        std::cout << "Blocks\t" << partition_config.k << std::endl;
        std::cout << "Seed\t" << partition_config.seed << std::endl;
        std::cout << "MLS stop global threshold\t" << partition_config.stop_mls_global_threshold << std::endl;
        std::cout << "MLS stop local threshold\t" << partition_config.stop_mls_local_threshold << std::endl;
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        // pin main thread to core
        parallel::PinToCore(partition_config.main_core);
        parallel::g_thread_pool.Resize(partition_config.num_threads - 1);

        std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
        if (partition_config.label_propagation_refinement) {
                std::cout << "Algorithm\t" << partition_config.configuration << std::endl;
                std::cout << "Block size\t" << partition_config.block_size << std::endl;
        } else {
                std::cout << "Algorithm\t" << partition_config.configuration + (partition_config.lp_before_local_search ? "_lp" : "") << std::endl;
        }
        std::cout << "Num threads\t" << partition_config.num_threads << std::endl;
        std::cout << "Max number of moves\t" << partition_config.max_number_of_moves << std::endl;

        if (partition_config.kway_all_boundary_nodes_refinement) {
                std::cout << "Refinement boundary\tall" << std::endl;
        } else {
                std::cout << "Refinement boundary\tpair" << std::endl;
        }

        if (partition_config.quotient_graph_two_way_refinement) {
                std::cout << "Two way refinement\ttrue" << std::endl;
        } else {
                std::cout << "Two way refinement\tfalse" << std::endl;
        }

        switch (partition_config.apply_move_strategy) {
                case ApplyMoveStrategy::LOCAL_SEARCH:
                        std::cout << "Move strategy\tlocal search" << std::endl;
                        break;
                case ApplyMoveStrategy::GAIN_RECALCULATION:
                        std::cout << "Move strategy\tgain recalculation" << std::endl;
                        break;
                case ApplyMoveStrategy::REACTIVE_VERTICES:
                        std::cout << "Move strategy\treactivate_vertices" << std::endl;
                        break;
                case ApplyMoveStrategy::SKIP:
                        std::cout << "Move strategy\tskip" << std::endl;
                        break;
        }

        switch (partition_config.kway_stop_rule) {
                case KWayStopRule::KWAY_SIMPLE_STOP_RULE:
                        std::cout << "Kway stop rule\tsimple" << std::endl;
                        break;
                case KWayStopRule::KWAY_ADAPTIVE_STOP_RULE:
                        std::cout << "Kway stop rule\tadaptive" << std::endl;
                        break;
                case KWayStopRule::KWAY_CHERNOFF_ADAPTIVE_STOP_RULE:
                        std::cout << "Kway stop rule\t" << kway_chernoff_adaptive_stop_rule::get_algo_name() << std::endl;
                        std::cout << "Stop probability\t" << partition_config.chernoff_stop_probability << std::endl;
                        std::cout << "Num gradient descent step\t" << partition_config.chernoff_gradient_descent_num_steps << std::endl;
                        std::cout << "Gradient descent step size\t" << partition_config.chernoff_gradient_descent_step_size << std::endl;
                        std::cout << "Min num step limit\t" << partition_config.chernoff_min_step_limit << std::endl;
                        std::cout << "Max num step limit\t" << partition_config.chernoff_max_step_limit << std::endl;
                        break;
        }

        // ***************************** perform partitioning ***************************************       
        t.restart();
        graph_partitioner partitioner;

        std::cout <<  "performing partitioning!"  << std::endl;
        if(partition_config.time_limit == 0) {
                partitioner.perform_partitioning(partition_config, G);
        } else {
                PartitionID* map = new PartitionID[G.number_of_nodes()];
                EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();
                while(t.elapsed() < partition_config.time_limit) {
                        partition_config.graph_allready_partitioned = false;
                        partitioner.perform_partitioning(partition_config, G);
                        EdgeWeight cut = qm.edge_cut(G);
                        if(cut < best_cut) {
                                best_cut = cut;
                                forall_nodes(G, node) {
                                        map[node] = G.getPartitionIndex(node);
                                } endfor
                        }
                }

                forall_nodes(G, node) {
                        G.setPartitionIndex(node, map[node]);
                } endfor
        }

        if( partition_config.kaffpa_perfectly_balance ) {
                double epsilon                         = partition_config.imbalance/100.0;
                partition_config.upper_bound_partition = (1+epsilon)*ceil(partition_config.largest_graph_weight/(double)partition_config.k);

                complete_boundary boundary(&G);
                boundary.build();

                cycle_refinement cr;
                cr.perform_refinement(partition_config, G, boundary);
        }
        // ******************************* done partitioning *****************************************       
        ofs.close();
        std::cout.rdbuf(backup);
        std::cout << "imbalance\t" << partition_config.imbalance / 100.0 << std::endl;
        std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;
       
        // output some information about the partition that we have computed
        Gain cut = qm.edge_cut(G);
        std::cout << "cut \t\t"         << cut                            << std::endl;
        if (partition_config.input_partition != "") {
                std::cout << "input partition cut\t" << input_partition_cut << std::endl;
                std::cout << "improvement\t" << input_partition_cut - cut << std::endl;
        }
        std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
        std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
        std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
        std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;

        if (!partition_config.label_propagation_refinement) {
                std::cout << "Two way refinement:" << std::endl;
                quotient_graph_refinement::print_full_statistics();
                std::cout << std::endl;

                std::cout << "Local search statistics:" << std::endl;
                if (partition_config.parallel_multitry_kway) {
                        parallel::multitry_kway_fm::print_full_statistics();
                        if (partition_config.input_partition != "") {
                                ALWAYS_ASSERT(
                                        parallel::multitry_kway_fm::get_performed_gain() == input_partition_cut - cut);
                        }
                } else {
                        multitry_kway_fm::print_full_statistics();
                }
                std::cout << std::endl;
        }

#ifdef PERFORMANCE_STATISTICS
        parallel::thread_data_refinement_core::nodes_partitions_hash_table::print_statistics();
#endif

#ifdef OUTPUT_GLOBAL_STAT
        kway_adaptive_stop_rule::dump_statistics();
        kway_chernoff_adaptive_stop_rule::dump_statistics();

        uint32_t ca_g_better = 0;
        int ca_g_better_val = 0;

        uint32_t a_g_better = 0;
        int a_g_better_val = 0;

        uint32_t ca_m_better = 0;
        uint32_t ca_m_better_val = 0;

        uint32_t a_m_better = 0;
        uint32_t a_m_better_val = 0;

        uint32_t equal = 0;

        for (size_t i = 0; i < kway_adaptive_stop_rule::m_stat_movements.size(); ++i) {
                uint32_t ca_m = kway_chernoff_adaptive_stop_rule::m_stat_movements[i];
                int ca_g = kway_chernoff_adaptive_stop_rule::m_stat_gains[i];

                uint32_t a_m = kway_adaptive_stop_rule::m_stat_movements[i];
                int a_g = kway_adaptive_stop_rule::m_stat_gains[i];

                if (a_g > ca_g) {
                    ++a_g_better;
                    a_g_better_val += a_g - ca_g;
                }
                if (a_g < ca_g) {
                    ++ca_g_better;
                    ca_g_better_val += ca_g - a_g;
                }
                if (a_g == ca_g) {
                    if (a_m < ca_m) {
                        ++a_m_better;
                        a_m_better_val += ca_m - a_m;
                    }
                    if (a_m > ca_m) {
                        ++ca_m_better;
                        ca_m_better_val += a_m - ca_m;
                    }
                    if (a_m == ca_m)
                        ++equal;
                }
        }

        std::cout << "Chernoff better by gain\t" << ca_g_better << "\ton\t" << ca_g_better_val << std::endl;
        std::cout << "Adaptive better by gain\t" << a_g_better << "\ton\t" << a_g_better_val << std::endl;
        std::cout << "Chernoff better by movement\t" << ca_m_better << "\ton\t" << ca_m_better_val << std::endl;
        std::cout << "Adaptive better by movement\t" << a_m_better << "\ton\t" << a_m_better_val << std::endl;
        std::cout << "Equal\t" << equal << std::endl;
#endif
        // write the partition to the disc
        if (!partition_config.filename_output.empty()) {
                graph_io::writePartition(G, partition_config.filename_output);
        }

        return 0;
}
