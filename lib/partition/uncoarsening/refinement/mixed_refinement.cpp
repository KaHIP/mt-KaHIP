/******************************************************************************
 * mixed_refinement.cpp 
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

#include "cycle_improvements/cycle_refinement.h"
#include "kway_graph_refinement/kway_graph_refinement.h"
#include "kway_graph_refinement/multitry_kway_fm.h"
#include "label_propagation_refinement/label_propagation_refinement.h"
#include "mixed_refinement.h"
#include "quotient_graph_refinement/quotient_graph_refinement.h"

#include "data_structure/parallel/time.h"

mixed_refinement::mixed_refinement() {
        multitry_kway_fm::reset_statistics();
}

mixed_refinement::~mixed_refinement() {

}

EdgeWeight mixed_refinement::perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary) {
        refinement* refine              = new quotient_graph_refinement();
        refinement* kway                = new kway_graph_refinement();
        cycle_refinement* cycle_refine  = new cycle_refinement();

        EdgeWeight overall_improvement = 0; 
        //call refinement
        if(config.no_change_convergence) {
                bool sth_changed = true;
                while(sth_changed) {
                        EdgeWeight improvement = 0;
                        if(config.corner_refinement_enabled) {
                                improvement += kway->perform_refinement(config, G, boundary);
                        }

                        if(!config.quotient_graph_refinement_disabled) {
                                improvement += refine->perform_refinement(config, G, boundary);
                        }

                        overall_improvement += improvement;
                        sth_changed = improvement != 0;
                }

        } else {
                if(config.corner_refinement_enabled && !config.parallel_multitry_kway) {
                        CLOCK_START;
                        quality_metrics qm;
                        EdgeWeight old_cut = 0;
                        if (config.check_cut) {
                                old_cut = qm.edge_cut(G);
                                std::cout << "before\t" << old_cut << std::endl;
                        }
                        EdgeWeight improvement = kway->perform_refinement(config, G, boundary);
                        overall_improvement += improvement;

                        if (config.check_cut) {
                                EdgeWeight new_cut = qm.edge_cut(G);
                                std::cout << "after\t" << new_cut << std::endl;
                                std::cout << "improvement\t" << improvement << std::endl;
                                ALWAYS_ASSERT(old_cut - new_cut == improvement);
                        }

                        CLOCK_END("Kway refinement");
                }

//                if(config.fastmultitry) {
//                        CLOCK_START;
//                        std::unique_ptr<multitry_kway_fm> multitry_kway = get_multitry_kway_fm_instance(config, G,
//                                                                                                        boundary);
//                        quality_metrics qm;
//                        EdgeWeight old_cut;
//                        double old_balance;
//                        if (config.check_cut) {
//                                old_cut = qm.edge_cut(G);
//                                old_balance = qm.balance(G);
//                                std::cout << "before cut\t" << old_cut << std::endl;
//                                std::cout << "balance before\t" << old_balance << std::endl;
//                        }
//
//
//                        EdgeWeight improvement = multitry_kway->perform_refinement(config, G, boundary,
//                                                                                 config.global_multitry_rounds, true,
//                                                                                 config.kway_adaptive_limits_alpha);
//                        overall_improvement += improvement;
//
//                        if (config.check_cut) {
//                                EdgeWeight new_cut = qm.edge_cut(G);
//                                double new_balance = qm.balance(G);
//                                std::cout << "after cut\t" << new_cut << std::endl;
//                                std::cout << "after balance\t" << new_balance << std::endl;
//                                std::cout << "improvement\t" << improvement << std::endl;
//                                ALWAYS_ASSERT(old_cut - new_cut == improvement);
//                        }
//
//                        CLOCK_END("Multitry kway refinement");
//                }

                if(!config.quotient_graph_refinement_disabled) {
                        CLOCK_START;
//                        if (!config.kway_all_boundary_nodes_refinement) {
//                                overall_improvement += refine->perform_refinement(config, G, boundary);
//                        } else {
//                                overall_improvement += refine->perform_refinement_all(config, G, boundary);
//                        }
                        overall_improvement += refine->perform_refinement(config, G, boundary);
                        CLOCK_END("Quotient graph refinement");
                }

                if(config.kaffpa_perfectly_balanced_refinement) {
                        CLOCK_START;
                        overall_improvement += cycle_refine->perform_refinement(config, G, boundary);
                        CLOCK_END("Cycle refinement");
                }
        }

        delete refine;
        delete kway;
        delete cycle_refine;

        return overall_improvement;
}

