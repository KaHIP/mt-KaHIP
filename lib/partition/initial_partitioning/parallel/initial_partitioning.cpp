#include "data_structure/parallel/random.h"
#include "data_structure/parallel/thread_pool.h"
#include "data_structure/parallel/time.h"
#include "graph_partition_assertions.h"
#include "initial_partitioning/initial_partition_bipartition.h"
#include "initial_partitioning/parallel/initial_partitioning.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "tools/macros_assertions.h"
#include "timer.h"

#include <fstream>

namespace parallel {

initial_partitioning::initial_partitioning() {
}

initial_partitioning::~initial_partitioning() {
}

void initial_partitioning::perform_initial_partitioning(const PartitionConfig& config, graph_access& G) {
        ALWAYS_ASSERT(config.initial_partitioning_type == INITIAL_PARTITIONING_RECPARTITION);

        quality_metrics qm;
        std::cout << "Parallel initial partitioning" << std::endl;
        timer t;
        t.restart();

        uint32_t reps_to_do = (unsigned) std::max((int)ceil(config.initial_partitioning_repetitions/(double)log2(config.k)),2);
        reps_to_do = std::max(reps_to_do, config.num_threads);

        std::atomic<uint32_t> reps_done(0);

        auto task_impl = [&qm, &config, &reps_done, &reps_to_do] (graph_access& G, uint32_t id)
                -> std::pair<EdgeWeight, std::unique_ptr<int[]>>{
                initial_partition_bipartition partition;

                EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();
                std::unique_ptr<int[]> best_map = std::make_unique<int[]>(G.number_of_nodes());
                ALWAYS_ASSERT(!(config.graph_allready_partitioned && !config.omit_given_partitioning));

                std::unique_ptr<int[]> partition_map = std::make_unique<int[]>(G.number_of_nodes());

                parallel::random rnd(id + config.seed);
                if (!((config.graph_allready_partitioned && config.no_new_initial_partitioning) ||
                      config.omit_given_partitioning)) {
                        uint32_t rep = 0;
                        while ((rep = reps_done.fetch_add(1, std::memory_order_release) < reps_to_do)) {
                                uint32_t seed = rnd.random_number(0u, std::numeric_limits<uint32_t>::max());
                                PartitionConfig working_config = config;
                                working_config.combine = false;

                                partition.initial_partition(working_config, seed, G, partition_map.get());

                                EdgeWeight cur_cut = qm.edge_cut(G, partition_map.get());
                                if (cur_cut < best_cut) {
                                        best_map.swap(partition_map);
                                        best_cut = cur_cut;

                                        if (best_cut == 0) {
                                                break;
                                        }
                                }
                        }
                }
                return std::make_pair(best_cut, std::move(best_map));
        };

        auto task = [&] (uint32_t id) {
                if (id > 0) {
                        random_functions::setSeed(id + config.seed);
                        graph_access my_graph;
                        G.copy(my_graph);
                        return task_impl(my_graph, id);
                } else {
                        return task_impl(G, id);
                }
        };

        std::vector<std::future<std::pair<EdgeWeight, std::unique_ptr<int[]>>>> futures;
        futures.reserve(g_thread_pool.NumThreads());

        // start initial
        std::ofstream ofs;
        std::streambuf* backup = std::cout.rdbuf();
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf());

        for (uint32_t id = 0; id < g_thread_pool.NumThreads(); ++id) {
                futures.push_back(parallel::g_thread_pool.Submit(id, task, id + 1));
        }

        EdgeWeight best_cut;
        std::unique_ptr<int[]> best_map;

        std::tie(best_cut, best_map) = task(0);

        parallel::random rnd(config.seed);
        std::vector<std::pair<EdgeWeight, std::unique_ptr<int[]>>> cuts;
        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                cuts.push_back(future.get());
        });
        ofs.close();
        std::cout.rdbuf(backup);

        for (auto& cut : cuts) {
                EdgeWeight cur_cut = cut.first;
                std::unique_ptr<int[]> cur_map = std::move(cut.second);

                if (cur_cut < best_cut || (cur_cut == best_cut && rnd.bit())) {
                        PRINT(std::cout << "log>"
                                        << "improved the current initial partitiong from "
                                        << best_cut
                                        << " to " << cur_cut << std::endl;)
                        best_cut = cur_cut;
                        best_map = std::move(cur_map);
                }
        }

        G.set_partition_count(config.k);
        forall_nodes(G, n) {
                G.setPartitionIndex(n, best_map[n]);
        } endfor

        PRINT(std::cout << "initial partitioning took " << t.elapsed()                << std::endl;)
        PRINT(std::cout << "log>"                       << "current initial balance " << qm.balance(G) << std::endl;)

        ALWAYS_ASSERT(!config.initial_partition_optimize && !config.combine);

        if(!(config.graph_allready_partitioned && config.no_new_initial_partitioning)) {
                PRINT(std::cout << "finalinitialcut " << best_cut                         << std::endl;)
                PRINT(std::cout << "log>"             << "final current initial balance " << qm.balance(G) << std::endl;)
        }
        std::cout << "initial cut\t" << best_cut << std::endl;

        ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, G));
}

}