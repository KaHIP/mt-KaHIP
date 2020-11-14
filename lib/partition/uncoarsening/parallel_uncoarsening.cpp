#include "data_structure/parallel/time.h"
#include "partition/uncoarsening/parallel_uncoarsening.h"
#include "partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "partition/uncoarsening/refinement/parallel_kway_graph_refinement/multitry_kway_fm.h"
#include "tools/graph_partition_assertions.h"
#include "tools/quality_metrics.h"

#include <memory>

namespace parallel {

int uncoarsening::perform_uncoarsening_cut(const PartitionConfig& config, graph_hierarchy& hierarchy) {
        PartitionConfig cfg = config;
        std::unique_ptr<graph_access> coarsest(hierarchy.get_coarsest());
        PRINT(std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;)

        if (config.lp_before_local_search) {
                CLOCK_START;
                perform_label_propagation(cfg, *coarsest);
                CLOCK_END(">> Uncoarsening: Label propagation");
        }

        double factor = config.balance_factor;
        cfg.upper_bound_partition = ((!hierarchy.isEmpty()) * factor + 1.0) * config.upper_bound_partition;

        EdgeWeight improvement = 0;
        if (config.parallel_multitry_kway) {
                CLOCK_START;
                boundary_type boundary(*coarsest, cfg);
                boundary.construct_boundary();
                if (config.check_cut) {
                        boundary.check_boundary();
                }
                CLOCK_END(">> Build boundary");

                CLOCK_START_N;
                improvement += perform_multitry_kway(cfg, *coarsest, boundary);
                CLOCK_END(">> Refinement");
        }

        uint32_t hierarchy_deepth = hierarchy.size();
        std::vector<std::unique_ptr<graph_access>> graphs_to_delete;

        while (!hierarchy.isEmpty()) {
                CLOCK_START;
                graph_access* G = hierarchy.parallel_pop_finer_and_project();
                CLOCK_END("Projection");

                PRINT(std::cout << "log>" << "unrolling graph with " << G->number_of_nodes() << std::endl;)

                if (config.lp_before_local_search) {
                        CLOCK_START_N;
                        perform_label_propagation(cfg, *G);
                        CLOCK_END(">> Uncoarsening: Label propagation");
                }

                //call refinement
                double cur_factor = factor / (hierarchy_deepth - hierarchy.size());
                cfg.upper_bound_partition = ((!hierarchy.isEmpty()) * cur_factor + 1.0) * config.upper_bound_partition;
                PRINT(std::cout << "cfg upperbound " << cfg.upper_bound_partition << std::endl;)

                if (config.parallel_multitry_kway) {
                        CLOCK_START_N;
                        boundary_type boundary(*G, cfg);
                        boundary.construct_boundary();
                        if (config.check_cut) {
                                boundary.check_boundary();
                        }
                        CLOCK_END(">> Build boundary");

                        CLOCK_START_N;
                        improvement += perform_multitry_kway(cfg, *G, boundary);
                        CLOCK_END(">> Refinement");
                }

                ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, *G));

                if (!hierarchy.isEmpty()) {
                        graphs_to_delete.emplace_back(G);
                }
        }

        return improvement;
}

void uncoarsening::perform_label_propagation(PartitionConfig& config, graph_access& G) {
        quality_metrics qm;
        EdgeWeight old_cut = 0;
        if (config.check_cut) {
                old_cut = qm.edge_cut(G);
                std::cout << "before\t" << old_cut << std::endl;
                std::cout << "upper_bound_partition\t" << config.upper_bound_partition << std::endl;
                std::cout << "before balance\t" << qm.balance(G) << std::endl;
        }

        EdgeWeight changed = label_propagation_refinement().perform_refinement(config, G);

        if (config.check_cut) {
                EdgeWeight new_cut = qm.edge_cut(G);
                std::cout << "after\t" << new_cut << std::endl;
                std::cout << "upper_bound_partition\t" << config.upper_bound_partition << std::endl;
                std::cout << "after balance\t" << qm.balance(G) << std::endl;
                std::cout << "changed\t" << changed << std::endl;
        }
}

EdgeWeight uncoarsening::perform_multitry_kway(PartitionConfig& config, graph_access& G, boundary_type& boundary) {
        CLOCK_START;
        quality_metrics qm;
        EdgeWeight old_cut = 0;
        double old_balance = 0.0;
        if (config.check_cut) {
                old_cut = qm.edge_cut(G);
                old_balance = qm.balance(G);
                std::cout << "before cut\t" << old_cut << std::endl;
                std::cout << "balance before\t" << old_balance << std::endl;
        }

        parallel::multitry_kway_fm multitry_kway(config, G, boundary);
        EdgeWeight improvement = multitry_kway.perform_refinement(config, G, boundary, config.global_multitry_rounds,
                                                                  true, config.kway_adaptive_limits_alpha);

        if (config.check_cut) {
                EdgeWeight new_cut = qm.edge_cut(G);
                double new_balance = qm.balance(G);
                std::cout << "after cut\t" << new_cut << std::endl;
                std::cout << "after balance\t" << new_balance << std::endl;
                std::cout << "improvement\t" << improvement << std::endl;
                ALWAYS_ASSERT(old_cut - new_cut == improvement);
        }
        CLOCK_END("Multitry kway refinement");
        return improvement;
}

}
