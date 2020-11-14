/******************************************************************************
 * coarsening.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include <limits>
#include <sstream>

#include "coarsening.h"
#include "coarsening_configurator.h"
#include "contraction.h"
#include "data_structure/graph_hierarchy.h"
#include "definitions.h"
#include "edge_rating/edge_ratings.h"
#include "graph_io.h"
#include "matching/gpa/gpa_matching.h"
#include "matching/random_matching.h"
#include "data_structure/parallel/time.h"
#include "min_hash/hash_common_neighborhood.h"

coarsening::coarsening() {

}

coarsening::~coarsening() {

}

void coarsening::perform_coarsening(const PartitionConfig & partition_config, graph_access & G, graph_hierarchy & hierarchy) {

        NodeID no_of_coarser_vertices = G.number_of_nodes();
        NodeID no_of_finer_vertices   = G.number_of_nodes();

        edge_ratings rating(partition_config);
        CoarseMapping* coarse_mapping = NULL;

        graph_access* finer                      = &G;
        matching* edge_matcher                   = NULL;
        contraction* contracter                  = new contraction();
        PartitionConfig copy_of_partition_config = partition_config;

        std::unique_ptr<stop_rule> coarsening_stop_rule = get_stop_rule(G, copy_of_partition_config);

        coarsening_configurator coarsening_config;

        unsigned int level    = 0;
        bool contraction_stop = false;
        bool common_neighborhood_clustering = partition_config.common_neighborhood_clustering;

        do {
                graph_access* coarser = new graph_access();
                coarse_mapping        = new CoarseMapping();
                Matching edge_matching;
                NodePermutationMap permutation;

                coarsening_config.configure_coarsening(copy_of_partition_config, &edge_matcher, level);

                CLOCK_START;
                if (partition_config.matching_type != CLUSTER_COARSENING) {
                        rating.rate(*finer, level);
                }
                CLOCK_END(">> Rate");

                CLOCK_START_N;
                edge_matcher->match(copy_of_partition_config, *finer, edge_matching,
                                    *coarse_mapping, no_of_coarser_vertices, permutation);
                delete edge_matcher;
                CLOCK_END(">> Match or clustering");

                if (common_neighborhood_clustering) {
                        std::cout << "Number of coarse vertices before min_hash = " << no_of_coarser_vertices << std::endl;
                        CLOCK_START;
                        hash_common_neighborhood().match(copy_of_partition_config, *finer, edge_matching,
                                                         *coarse_mapping, no_of_coarser_vertices, permutation);

                        std::cout << ">> common_neighborhood_clustering: finer vertices = " << finer->number_of_nodes() << ", # of coarsed vertices = " << no_of_coarser_vertices << std::endl;
                        if (finer->number_of_nodes() == no_of_coarser_vertices) {
                                common_neighborhood_clustering = false;
                        }
                        CLOCK_END(">> Hashing neighbors");
                }

//                if (!copy_of_partition_config.accept_small_coarser_graphs && no_of_coarser_vertices < copy_of_partition_config.k * 1000) {
//                        std::cout << "Do not accept this clustering. The number of vertices " << no_of_coarser_vertices << " < k * " << 1000 << std::endl;
//                        std::cout << "Number of vertices = " << no_of_coarser_vertices << std::endl;
//                        break;
//                }

                CLOCK_START_N;
                if(partition_config.graph_allready_partitioned) {
                        contracter->contract_partitioned(copy_of_partition_config, *finer, *coarser, edge_matching, 
                                                         *coarse_mapping, no_of_coarser_vertices, permutation);
                } else {
                        contracter->contract(copy_of_partition_config, *finer, *coarser, edge_matching,
                                             *coarse_mapping, no_of_coarser_vertices, permutation);
                }
                CLOCK_END(">> Contract");

                hierarchy.push_back(finer, coarse_mapping);

                contraction_stop = coarsening_stop_rule->stop(no_of_finer_vertices, *coarser);

                no_of_finer_vertices = no_of_coarser_vertices;
                std::cout <<  "no of coarser vertices " << coarser->number_of_nodes() <<  " and no of edges " <<  coarser->number_of_edges() << std::endl;

                finer = coarser;

                level++;

        } while (contraction_stop);

        hierarchy.push_back(finer, NULL); // append the last created level

        delete contracter;
}


