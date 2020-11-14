/******************************************************************************
 * multitry_kway_fm.h 
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

#ifndef MULTITRY_KWAYFM_PVGY97EW
#define MULTITRY_KWAYFM_PVGY97EW

#include  <vector>

#include "data_structure/parallel/algorithm.h"
#include "definitions.h"
#include "kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/refinement.h"

#include <memory>


class multitry_kway_fm {
public:
        multitry_kway_fm();

        virtual ~multitry_kway_fm();

        virtual int perform_refinement(PartitionConfig& config, graph_access& G,
                               complete_boundary& boundary, unsigned rounds,
                               bool init_neighbors, unsigned alpha);

        virtual int perform_refinement_around_parts(PartitionConfig& config, graph_access& G,
                                            complete_boundary& boundary, bool init_neighbors,
                                            unsigned alpha,
                                            PartitionID& lhs, PartitionID& rhs,
                                            std::unordered_map <PartitionID, PartitionID>& touched_blocks);

        static void reset_statistics() {
//                time_setup_start_nodes = 0.0;
//                time_local_search = 0.0;
//                time_generate_moves = 0.0;
//                tried_movements = 0;
//                kway_graph_refinement_commons::scaned_movements = 0;
//                kway_graph_refinement_commons::num_part_accesses = 0;
        }

        static void print_full_statistics() {
//                std::cout << "Time setup start nodes\t" << time_setup_start_nodes << " s" << std::endl;
//                std::cout << "Time local search\t" << time_local_search << " s" << std::endl;
//                std::cout << "Time generate moves\t" << time_generate_moves << " s" << std::endl;
//                std::cout << "Total tried moves\t" << tried_movements << std::endl;
//                std::cout << "Total scanned neighbours\t" << kway_graph_refinement_commons::scaned_movements << std::endl;
//                std::cout << "Number of partition accesses\t" << kway_graph_refinement_commons::num_part_accesses << std::endl;
        }

//        static double time_setup_start_nodes;
//        static double time_local_search;
//        static double time_generate_moves;
//        static uint32_t tried_movements;
private:

        int start_more_locallized_search(PartitionConfig& config, graph_access& G,
                                         complete_boundary& boundary,
                                         bool init_neighbors,
                                         bool compute_touched_blocks,
                                         std::unordered_map<PartitionID, PartitionID>& touched_blocks,
                                         std::vector<NodeID>& todolist);

        kway_graph_refinement_commons* commons;
};

std::unique_ptr<multitry_kway_fm> get_multitry_kway_fm_instance(PartitionConfig& config,
                                                                graph_access& G, complete_boundary& boundary);

#endif /* end of include guard: MULTITRY_KWAYFM_PVGY97EW  */
