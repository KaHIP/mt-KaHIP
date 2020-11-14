/******************************************************************************
 * quotient_graph_refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef QUOTIENT_GRAPH_REFINEMENT_A0Y1Y6LL
#define QUOTIENT_GRAPH_REFINEMENT_A0Y1Y6LL

#include "definitions.h"
#include "uncoarsening/refinement/refinement.h"

class quotient_graph_refinement : public refinement {
        public:
                quotient_graph_refinement( );
                virtual ~quotient_graph_refinement();

                EdgeWeight perform_refinement(PartitionConfig & config, graph_access & G,
                                              complete_boundary & boundary) override;

                EdgeWeight perform_refinement_all(PartitionConfig& config, graph_access& G,
                                                  complete_boundary& boundary) override;

                void setup_start_nodes(graph_access & G, 
                                       PartitionID partition, 
                                       boundary_pair & bp, 
                                       complete_boundary & boundary,  
                                       boundary_starting_nodes & start_nodes);

                static void print_full_statistics() {
                        //std::cout << "Total time two way\t" << total_time_two_way << std::endl;
                }
        private:
                //static double total_time_two_way;

                EdgeWeight perform_a_two_way_refinement(PartitionConfig & config, 
                                                        graph_access & G,
                                                        complete_boundary & boundary, 
                                                        boundary_pair & bp,
                                                        PartitionID & lhs, 
                                                        PartitionID & rhs,
                                                        NodeWeight & lhs_part_weight,
                                                        NodeWeight & rhs_part_weight,
                                                        EdgeWeight & cut,
                                                        bool & something_changed); 

};


#endif /* end of include guard: QUOTIENT_GRAPH_REFINEMENT_A0Y1Y6LL */
