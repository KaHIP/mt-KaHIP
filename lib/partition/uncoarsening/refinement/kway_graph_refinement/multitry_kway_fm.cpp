/******************************************************************************
 * multitry_kway_fm.cpp 
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

#include <algorithm>
#include <unordered_map>

#include "data_structure/parallel/time.h"
#include "kway_graph_refinement_core.h"
#include "kway_stop_rule.h"
#include "multitry_kway_fm.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"
#include "partition/uncoarsening/refinement/parallel_kway_graph_refinement/multitry_kway_fm.h"

//double multitry_kway_fm::time_setup_start_nodes(0);
//double multitry_kway_fm::time_local_search(0);
//double multitry_kway_fm::time_generate_moves(0);
//uint32_t multitry_kway_fm::tried_movements(0);

std::unique_ptr<multitry_kway_fm> get_multitry_kway_fm_instance(PartitionConfig& config,
                                                                graph_access& G, complete_boundary& boundary) {
//        if (config.parallel_multitry_kway) {
//                return std::make_unique<parallel::multitry_kway_fm>(config, G, boundary);
//        } else {
                return std::make_unique<multitry_kway_fm>();
//        }
}

multitry_kway_fm::multitry_kway_fm() {
}

multitry_kway_fm::~multitry_kway_fm() {

}

int multitry_kway_fm::perform_refinement(PartitionConfig& config, graph_access& G,
                                             complete_boundary& boundary, unsigned rounds,
                                             bool init_neighbors, unsigned alpha) {

        commons = kway_graph_refinement_commons::getInstance(config);

        unsigned tmp_alpha = config.kway_adaptive_limits_alpha;
        KWayStopRule tmp_stop = config.kway_stop_rule;
        config.kway_adaptive_limits_alpha = alpha;
        config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;

        int overall_improvement = 0;
        //while (true) {
        for (unsigned i = 0; i < rounds; i++) {
                CLOCK_START;
                boundary_starting_nodes start_nodes;
                boundary.setup_start_nodes_all(G, start_nodes);
                //time_setup_start_nodes += CLOCK_END_TIME;

                if (start_nodes.size() == 0) {
                        return 0;
                }// nothing to refine

                //now we do something with the start nodes
                //convert it into a list
                std::vector <NodeID> todolist;
                for (unsigned i = 0; i < start_nodes.size(); i++) {
                        todolist.push_back(start_nodes[i]);
                }

                CLOCK_START_N;
                std::unordered_map<PartitionID, PartitionID> touched_blocks;
                EdgeWeight improvement = start_more_locallized_search(config, G, boundary, init_neighbors, false,
                                                                      touched_blocks, todolist);
                //time_local_search += CLOCK_END_TIME;

                if (improvement == 0) {
                        break;
                }
                overall_improvement += improvement;

        }
        ASSERT_TRUE(overall_improvement >= 0);

        config.kway_adaptive_limits_alpha = tmp_alpha;
        config.kway_stop_rule = tmp_stop;

        return (int) overall_improvement;

}

int multitry_kway_fm::perform_refinement_around_parts(PartitionConfig& config, graph_access& G,
                                                          complete_boundary& boundary, bool init_neighbors,
                                                          unsigned alpha,
                                                          PartitionID& lhs, PartitionID& rhs,
                                                          std::unordered_map <PartitionID, PartitionID>& touched_blocks) {
        commons = kway_graph_refinement_commons::getInstance(config);

        unsigned tmp_alpha = config.kway_adaptive_limits_alpha;
        KWayStopRule tmp_stop = config.kway_stop_rule;
        config.kway_adaptive_limits_alpha = alpha;
        //config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;
        int overall_improvement = 0;

        for (unsigned i = 0; i < config.local_multitry_rounds; i++) {
                CLOCK_START;
                boundary_starting_nodes start_nodes;
                boundary.setup_start_nodes_around_blocks(G, lhs, rhs, start_nodes);
                //time_setup_start_nodes += CLOCK_END_TIME;

                if (start_nodes.size() == 0) {
                        return 0;
                }// nothing to refine

                //now we do something with the start nodes
                std::vector <NodeID> todolist;
                for (unsigned i = 0; i < start_nodes.size(); i++) {
                        todolist.push_back(start_nodes[i]);
                }

                CLOCK_START_N;
                EdgeWeight improvement = start_more_locallized_search(config, G, boundary,
                                                                      init_neighbors, true,
                                                                      touched_blocks, todolist);
                //time_local_search += CLOCK_END_TIME;

                if (improvement == 0) break;

                overall_improvement += improvement;
        }

        config.kway_adaptive_limits_alpha = tmp_alpha;
        config.kway_stop_rule = tmp_stop;
        ASSERT_TRUE(overall_improvement >= 0);
        return (int) overall_improvement;
}

int multitry_kway_fm::start_more_locallized_search(PartitionConfig& config, graph_access& G,
                                                   complete_boundary& boundary, bool init_neighbors,
                                                   bool compute_touched_blocks,
                                                   std::unordered_map <PartitionID, PartitionID>& touched_blocks,
                                                   std::vector <NodeID>& todolist) {
#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
        std::sort(todolist.begin(), todolist.end());
#endif

        random_functions::permutate_vector_good(todolist, false);
        commons = kway_graph_refinement_commons::getInstance(config);

        kway_graph_refinement_core refinement_core;
        int local_step_limit = 0;

        vertex_moved_hashtable moved_idx;
        unsigned idx = todolist.size() - 1;
        int overall_improvement = 0;

        CLOCK_START;
        while (!todolist.empty()) {
                int random_idx = random_functions::nextInt(0, idx);
                NodeID node = todolist[random_idx];

                PartitionID maxgainer;
                EdgeWeight extdeg = 0;
                commons->compute_gain(G, node, maxgainer, extdeg);

                if (moved_idx.find(node) == moved_idx.end() && extdeg > 0) {
                        boundary_starting_nodes real_start_nodes;
                        real_start_nodes.push_back(node);

                        if (init_neighbors) {
                                forall_out_edges(G, e, node)
                                {
                                        //++kway_graph_refinement_commons::scaned_movements;
                                        NodeID target = G.getEdgeTarget(e);
                                        if (moved_idx.find(target) == moved_idx.end()) {
                                                extdeg = 0;
                                                commons->compute_gain(G, target, maxgainer, extdeg);
                                                if (extdeg > 0) {
                                                        real_start_nodes.push_back(target);
                                                }
                                        }
                                }
                                endfor
                        }
                        int improvement = 0;
                        uint32_t tried_movements = 0;
                        if (compute_touched_blocks) {
                                std::tie(improvement, tried_movements) = refinement_core.single_kway_refinement_round(config, G,
                                                                                           boundary, real_start_nodes,
                                                                                           local_step_limit, moved_idx,
                                                                                           touched_blocks);
                                if (improvement < 0) {
                                        std::cout << "buf error improvement < 0" << std::endl;
                                }
                        } else {
                                std::tie(improvement, tried_movements) = refinement_core.single_kway_refinement_round(config, G,
                                                                                           boundary, real_start_nodes,
                                                                                           local_step_limit, moved_idx);
                                if (improvement < 0) {
                                        std::cout << "buf error improvement < 0" << std::endl;
                                }
                        }

                        overall_improvement += improvement;
                        //this->tried_movements += tried_movements;
                }

                if (!config.kway_all_boundary_nodes_refinement && moved_idx.size() > 0.05 * G.number_of_nodes()) break;
                std::swap(todolist[random_idx], todolist[idx--]);
                todolist.pop_back();
        }
        //auto t = CLOCK_END_TIME;
        //time_generate_moves += t;

        return overall_improvement;
}
