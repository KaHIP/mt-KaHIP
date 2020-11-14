/******************************************************************************
 * kway_graph_refinement_core.cpp
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
#include <tbb/concurrent_queue.h>

#include "data_structure/parallel/time.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "tools/random_functions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_stop_rule.h"
#include "uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_core.h"

namespace parallel {
constexpr unsigned int kway_graph_refinement_core::sentinel;
constexpr int kway_graph_refinement_core::signed_sentinel;

std::tuple<EdgeWeight, int, uint32_t>
kway_graph_refinement_core::single_kway_refinement_round(thread_data_refinement_core& td) {
        return single_kway_refinement_round_internal(td);
}

std::tuple<EdgeWeight, int, uint32_t>
kway_graph_refinement_core::single_kway_refinement_round_internal(thread_data_refinement_core& td) {
        auto& queue = td.queue;

        queue->clear();
        td.move_to->clear();

        init_queue_with_boundary(td, queue);

        if (queue->empty() || td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
                td.transpositions.push_back(sentinel);
                td.from_partitions.push_back(sentinel);
                td.to_partitions.push_back(sentinel);
                td.gains.push_back(signed_sentinel);

                return std::make_tuple(0, -1, 0);
        }

        int max_number_of_swaps = td.config.max_number_of_moves != -1 ? td.config.max_number_of_moves : td.G.number_of_nodes();

        EdgeWeight cut = std::numeric_limits<int>::max() / 2; // so we dont need to compute the edge cut
        EdgeWeight initial_cut = cut;

        //roll forwards
        EdgeWeight best_cut = cut;
        NodeID max_rnd = 0;
        int number_of_swaps = 0;
        uint32_t movements = 0;

        std::unique_ptr<kway_stop_rule> stopping_rule;

        switch (td.config.kway_stop_rule) {
                case KWAY_SIMPLE_STOP_RULE:
                        stopping_rule = std::make_unique<kway_simple_stop_rule>(td.config);
                        break;
                case KWAY_ADAPTIVE_STOP_RULE:
                        stopping_rule = std::make_unique<kway_adaptive_stop_rule>(td.config);
                        break;
                case KWAY_CHERNOFF_ADAPTIVE_STOP_RULE:
                        stopping_rule = std::make_unique<kway_chernoff_adaptive_stop_rule>(td.config);
                        break;
        }

        int previously_moved = (int) td.transpositions.size();

        // minus 1 for sentinel
        int min_cut_index = previously_moved - 1;
        for (number_of_swaps = 0, movements = 0; (int32_t) movements < max_number_of_swaps; movements++, number_of_swaps++) {
                if (queue->empty()) {
                        ++td.stop_empty_queue;
                        break;
                }

                uint32_t local_min_cut_index = (uint32_t) (min_cut_index - previously_moved >= 0 ? min_cut_index -
                                                                                                   previously_moved
                                                                                                : 0);
                if (stopping_rule->search_should_stop(local_min_cut_index, (uint32_t) number_of_swaps,
                                                      (uint32_t) td.step_limit)) {
                        ++td.stop_stopping_rule;
                        break;
                }

                Gain gain = queue->maxValue();
                NodeID node = queue->deleteMax();

                PartitionID from = td.get_local_partition(node);
                PartitionID to = td.get_local_partition_to_move(node);
                td.remove_local_partition_to_move(node);

                if (td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
                        break;
                }

                bool successfull = local_move_node(td, node, from, to, queue, gain);

                if (td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
                        break;
                }

                if (successfull) {
                        ++td.accepted_movements;
                        cut -= gain;
                        stopping_rule->push_statistics(gain);

#ifdef COMPARE_WITH_SEQUENTIAL_KAHIP
                        bool accept_equal = random_functions::nextBool();
#else
                        //bool accept_equal = td.rnd.bit();
                        NodeID cur_rnd = 0;
#endif
                        if (cut < best_cut || (cut == best_cut && max_rnd < (cur_rnd = td.rnd.random_number<NodeID>()) )) {
                                if (cut < best_cut) {
                                        cur_rnd = td.rnd.random_number<NodeID>();
                                        stopping_rule->reset_statistics();
                                }
                                best_cut = cut;
                                max_rnd = cur_rnd;
                                min_cut_index = previously_moved + number_of_swaps;
                        }
                        td.from_partitions.push_back(from);
                        td.to_partitions.push_back(to);
                        td.transpositions.push_back(node);
                        td.gains.push_back(gain);

                        ALWAYS_ASSERT(min_cut_index < (int64_t) td.transpositions.size());
                } else {
                        number_of_swaps--; //because it wasnt swaps
                }
        }

        if ((int32_t) movements == max_number_of_swaps) {
                ++td.stop_max_number_of_swaps;
        }

        uint32_t unrolled_moves = unroll_moves(td, min_cut_index);
        td.accepted_movements -= unrolled_moves;
        td.nodes_partitions->clear();

        td.transpositions.push_back(sentinel);
        td.from_partitions.push_back(sentinel);
        td.to_partitions.push_back(sentinel);
        td.gains.push_back(signed_sentinel);

        return std::make_tuple(initial_cut - best_cut, min_cut_index, movements);
}

uint32_t kway_graph_refinement_core::unroll_moves(thread_data_refinement_core& td, int min_cut_index) const {
        uint32_t unrolled_moves = 0;
        // unrolled_moves <  td.transpositions.size() - (min_cut_index + 1)
        while (unrolled_moves + min_cut_index + 1 < td.transpositions.size()) {
                size_t index = td.transpositions.size() - 1 - unrolled_moves;
                NodeID node = td.transpositions[index];
                PartitionID from = td.from_partitions[index];
                PartitionID to = td.to_partitions[index];
                local_move_back_node(td, node, from, to);

                ++unrolled_moves;
        }

        return unrolled_moves;
}

EdgeWeight kway_graph_refinement_core::apply_moves(uint32_t num_threads,
                                                   Cvector <thread_data_refinement_core>& threads_data,
                                                   std::vector<NodeID>& reactivated_vertices) const {

        EdgeWeight overall_gain = 0;
        //CLOCK_START;
        std::vector<std::future<void>> futures;
        std::atomic<uint32_t> thread_id(0);
        std::vector<AtomicWrapper<uint32_t>> offsets(num_threads, 0);

        if (parallel::g_thread_pool.NumThreads() > 0) {
                CLOCK_START;
                threads_data[0].get().boundary.begin_movements();
                threads_data[0].get().time_move_nodes_change_boundary += CLOCK_END_TIME;

                futures = prepare_boundary(num_threads, threads_data, thread_id, offsets);
        }
        //CLOCK_END("prepare boundary movements 1");

        for (size_t id = 0; id < num_threads; ++id) {
                overall_gain += apply_moves(threads_data[id].get(), reactivated_vertices, futures);
        }

        //CLOCK_START_N;
        {
                CLOCK_START;
                std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                        future.get();
                });
                threads_data[0].get().time_move_nodes_change_boundary += CLOCK_END_TIME;
        }
        //CLOCK_END("prepare boundary movements 2");
        return overall_gain;
}

std::vector<std::future<void>> kway_graph_refinement_core::prepare_boundary(uint32_t num_threads,
                                                                            Cvector<thread_data_refinement_core>& threads_data,
                                                                            std::atomic<uint32_t>& thread_id,
                                                                            std::vector<AtomicWrapper<uint32_t>>& offsets) const {

        std::vector<std::future<void>> futures;
        futures.reserve(parallel::g_thread_pool.NumThreads());
        auto task = [&, num_threads]() {
                while (true) {
                        uint32_t cur_thread_id = thread_id.load(std::memory_order_relaxed);

                        if (cur_thread_id >= num_threads) {
                                break;
                        }
                        auto& td = threads_data[cur_thread_id].get();

                        uint32_t i = offsets[cur_thread_id].fetch_add(1, std::memory_order_relaxed);
                        if (i >= td.min_cut_indices.size()) {
                                thread_id.compare_exchange_strong(cur_thread_id, cur_thread_id + 1,
                                                                  std::memory_order_release);
                                continue;
                        }

                        if (td.min_cut_indices[i].first == -1) {
                                continue;
                        }

                        uint32_t end = (uint32_t) td.min_cut_indices[i].first;

                        uint32_t begin = 0;
                        if (i > 0) {
                                begin = (uint32_t) td.min_cut_indices[i - 1].second + 1;
                        }

                        for (NodeID index = begin; index <= end; ++index) {
                                NodeID node = td.transpositions[index];

                                td.boundary.add_vertex_to_check(node);

                                if (m_maintain_complete_boundary) {
                                        forall_out_edges(td.G, e, node){
                                                NodeID target = td.G.getEdgeTarget(e);
                                                PartitionID part = td.G.getPartitionIndex(target);
                                                if (part == td.to_partitions[index] || part == td.from_partitions[index]) {
                                                        td.boundary.add_vertex_to_check(target);
                                                }
                                        } endfor
                                }
                        }
                }
        };

        {
                CLOCK_START;
                for (size_t i = 0; i < parallel::g_thread_pool.NumThreads(); ++i) {
                        futures.push_back(parallel::g_thread_pool.Submit(i, task));
                }
                threads_data[0].get().time_move_nodes_change_boundary += CLOCK_END_TIME;
        }

        return futures;
}

void kway_graph_refinement_core::update_boundary(thread_data_refinement_core& td) const {
        //CLOCK_START;
        if (parallel::g_thread_pool.NumThreads() > 0) {
                CLOCK_START;
                td.boundary.finish_movements();
                td.time_move_nodes_change_boundary += CLOCK_END_TIME;
        }
        //CLOCK_END("update boundary movements");
}

EdgeWeight kway_graph_refinement_core::apply_moves(thread_data_refinement_core& td,
                                                   std::vector<NodeID>& reactivated_vertices,
                                                   std::vector<std::future<void>>& futures) const {
        CLOCK_START;
        ALWAYS_ASSERT(td.transpositions.size() == td.from_partitions.size());
        ALWAYS_ASSERT(td.transpositions.size() == td.to_partitions.size());
        ALWAYS_ASSERT(td.transpositions.size() == td.gains.size());
        td.transpositions_size += td.transpositions.size();

        auto min_cut_iter = td.min_cut_indices.begin();
        EdgeWeight cut_improvement = 0;
        Gain total_expected_gain = 0;
        std::vector<NodeID> transpositions;
        std::vector<PartitionID> from_partitions;
        std::vector<EdgeWeight> gains;
        uint32_t aff = 0;

        for (int index = 0; index < (int) td.transpositions.size(); ++index) {
                int min_cut_index = min_cut_iter->first;
                int next_index = min_cut_iter->second;
                ++min_cut_iter;

                if (min_cut_index == -1) {
                        index = next_index;
                        continue;
                }

                int best_total_gain = 0;
                int total_gain = 0;
                NodeID max_rnd = 0;

                transpositions.reserve(min_cut_index - index + 1);
                from_partitions.reserve(min_cut_index - index + 1);
                gains.reserve(min_cut_index - index + 1);
                transpositions.clear();
                from_partitions.clear();
                gains.clear();

                while (index <= min_cut_index) {
                        NodeID node = td.transpositions[index];
                        PartitionID expected_from = td.from_partitions[index];
                        PartitionID expected_to = td.to_partitions[index];
                        EdgeWeight expected_gain = td.gains[index];

                        PartitionID from = td.G.getPartitionIndex(node);
                        PartitionID to;
                        EdgeWeight gain = td.compute_gain_actual(node, from, to, expected_to);

                        bool same_move = true;
                        if (expected_from != from || expected_to != to || expected_gain != gain) {
                                same_move = false;
                                ++td.affected_movements;
                                ++aff;
                                if (g_thread_pool.NumThreads() > 0 && m_maintain_complete_boundary) {
                                        CLOCK_START;
                                        auto task = [&td] (NodeID cur_node, PartitionID to, PartitionID from) {
                                                forall_out_edges(td.G, e, cur_node) {
                                                        NodeID target = td.G.getEdgeTarget(e);
                                                        PartitionID part = td.G.getPartitionIndex(target);
                                                        if (part == from || part == to) {
                                                                td.boundary.add_vertex_to_check(target);
                                                        }
                                                } endfor
                                        };
                                        uint32_t thread_id = td.rnd.random_number(0u, (uint32_t) parallel::g_thread_pool.NumThreads() - 1);
                                        futures.push_back(parallel::g_thread_pool.Submit(thread_id, task, node, to, from));
                                        td.time_move_nodes_change_boundary += CLOCK_END_TIME;
                                }
                        }

                        if (to == INVALID_PARTITION) {
                                ALWAYS_ASSERT(aff > 0);
                                ++index;
                                continue;
                        }

                        total_expected_gain += expected_gain;

                        bool success = relaxed_move_node(td, node, from, to);

                        if (success) {
                                transpositions.push_back(node);
                                from_partitions.push_back(from);
                                gains.push_back(gain);

                                cut_improvement += gain;
                                total_gain += gain;

                                NodeID cur_rnd = 0;
                                if (total_gain > best_total_gain || (total_gain == best_total_gain && ( (cur_rnd = td.rnd.random_number<NodeID>()) > max_rnd || same_move))) {
                                        if (total_gain > best_total_gain) {
                                                cur_rnd = td.rnd.random_number<NodeID>();
                                        }

                                        best_total_gain = total_gain;
                                        max_rnd = cur_rnd;

                                        for (size_t i = 0; i < transpositions.size(); ++i) {
                                                NodeID node = transpositions[i];
                                                reactivated_vertices.push_back(node);
                                        }

                                        from_partitions.clear();
                                        transpositions.clear();
                                        gains.clear();
                                }
                        }
                        ++index;
                }
                unroll_relaxed_moves(td, transpositions, from_partitions, gains, cut_improvement);

                index = next_index;
        }

        td.time_move_nodes += CLOCK_END_TIME;
        td.unperformed_gain += total_expected_gain - cut_improvement;
        td.performed_gain += cut_improvement;
        return cut_improvement;
}

void kway_graph_refinement_core::init_queue_with_boundary(thread_data_refinement_core& td,
                                                          std::unique_ptr<refinement_pq>& queue) {
        if (td.config.permutation_during_refinement == PERMUTATION_QUALITY_FAST) {
                random_functions::permutate_vector_fast(td.start_nodes, false);
        } else if (td.config.permutation_during_refinement == PERMUTATION_QUALITY_GOOD) {
                random_functions::permutate_vector_good(td.start_nodes, false);
        }

        for (NodeID node : td.start_nodes) {
                bool expected = false;
                if (td.moved_idx[node].compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                        PartitionID to;
                        EdgeWeight ext_degree;
                        //compute gain
                        PartitionID from = td.get_local_partition(node);

                        Gain gain = td.compute_gain(node, from, to, ext_degree);

                        if (td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
                                return;
                        }

                        queue->insert(node, gain);
                        td.set_local_partition_to_move(node, to);
                        td.moved.push_back(node);
                }
        }
}

inline bool kway_graph_refinement_core::relaxed_move_node(thread_data_refinement_core& td,
                                                          NodeID node,
                                                          PartitionID from,
                                                          PartitionID to) const {
        ALWAYS_ASSERT(td.G.getPartitionIndex(node) == from);

        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);

        if (td.boundary.get_block_weight(to) + this_nodes_weight >= td.config.upper_bound_partition) {
                return false;
        }

        if (td.boundary.get_block_size(from) == 1) {// assure that no block gets accidentally empty
                return false;
        }

        td.G.setPartitionIndex(node, to);

        if (parallel::g_thread_pool.NumThreads() == 0) {
                CLOCK_START;
                td.boundary.move(node, from, to);
                auto t = CLOCK_END_TIME;
                td.time_move_nodes_change_boundary += t;
        }

        td.boundary.set_block_size(from, td.boundary.get_block_size(from) - 1);
        td.boundary.set_block_size(to, td.boundary.get_block_size(to) + 1);
        td.boundary.set_block_weight(from, td.boundary.get_block_weight(from) - this_nodes_weight);
        td.boundary.set_block_weight(to, td.boundary.get_block_weight(to) + this_nodes_weight);

        return true;
}

void kway_graph_refinement_core::unroll_relaxed_moves(thread_data_refinement_core& td,
                                                      std::vector<NodeID>& transpositions,
                                                      std::vector<PartitionID>& from_partitions,
                                                      std::vector<Gain>& gains,
                                                      int& cut_improvement) const {
        size_t size = transpositions.size();
        for (size_t i = 0; i < transpositions.size(); ++i) {
                size_t index = size - i - 1;
                NodeID node = transpositions[index];
                PartitionID from = from_partitions[index];
                PartitionID to = td.G.getPartitionIndex(node);
                cut_improvement -= gains[index];
                relaxed_move_node_back(td, node, from, to);
        }
}

void kway_graph_refinement_core::relaxed_move_node_back(thread_data_refinement_core& td, NodeID node, PartitionID from,
                                                        PartitionID to) const {
        ALWAYS_ASSERT(td.G.getPartitionIndex(node) == to);
        td.G.setPartitionIndex(node, from);

        if (parallel::g_thread_pool.NumThreads() == 0) {
                CLOCK_START;
                td.boundary.move(node, to, from);
                auto t = CLOCK_END_TIME;
                td.time_move_nodes_change_boundary += t;
        }

        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);
        td.boundary.set_block_size(from, td.boundary.get_block_size(from) + 1);
        td.boundary.set_block_size(to, td.boundary.get_block_size(to) - 1);
        td.boundary.set_block_weight(from, td.boundary.get_block_weight(from) + this_nodes_weight);
        td.boundary.set_block_weight(to, td.boundary.get_block_weight(to) - this_nodes_weight);
}

inline bool kway_graph_refinement_core::local_move_back_node(thread_data_refinement_core& td,
                                                             NodeID node,
                                                             PartitionID from,
                                                             PartitionID to) const {
        // td.set_local_partition(node, from);
        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);

//        td.parts_weights[from].get().fetch_add(this_nodes_weight, std::memory_order_relaxed);
//        td.parts_weights[to].get().fetch_sub(this_nodes_weight, std::memory_order_relaxed);
//        td.parts_sizes[to].get().fetch_sub(1, std::memory_order_relaxed);
//        td.parts_sizes[from].get().fetch_add(1, std::memory_order_relaxed);

        td.parts_weights[from] += this_nodes_weight;
        td.parts_weights[to] -= this_nodes_weight;
        --td.parts_sizes[to];
        ++td.parts_sizes[from];


        return true;
};

inline bool kway_graph_refinement_core::local_move_node(thread_data_refinement_core& td,
                                                       NodeID node,
                                                       PartitionID from,
                                                       PartitionID to,
                                                       std::unique_ptr<refinement_pq>& queue, Gain gain) {

//        EdgeWeight node_ext_deg;
//        PartitionID expected_to;
//        Gain expected_gain = td.compute_gain(node, from, expected_to, node_ext_deg);

//        if (td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
//                return false;
//        }

//        ALWAYS_ASSERT(expected_gain == gain);
//        ALWAYS_ASSERT(to != INVALID_PARTITION);

        NodeWeight this_nodes_weight = td.G.getNodeWeight(node);

//        if (td.parts_sizes[from].get().load(std::memory_order_relaxed) == 1) {
//                return false;
//        }
//
//        NodeWeight part_weight = td.parts_weights[to].get().load(std::memory_order_relaxed);
//
//        do {
//                if (part_weight + this_nodes_weight >= td.config.upper_bound_partition) {
//                        return false;
//                }
//        } while (!td.parts_weights[to].get().compare_exchange_weak(part_weight,
//                                                                   part_weight + this_nodes_weight,
//                                                                   std::memory_order_relaxed));

        if (td.parts_sizes[from] == 1) {
                return false;
        }

        NodeWeight part_weight = td.parts_weights[to];

        if (part_weight + this_nodes_weight >= td.config.upper_bound_partition) {
                return false;
        }
        td.parts_weights[to] = part_weight + this_nodes_weight;


        td.set_local_partition(node, to);
//        td.parts_weights[from].get().fetch_sub(this_nodes_weight, std::memory_order_relaxed);
//        td.parts_sizes[to].get().fetch_add(1, std::memory_order_relaxed);
//        td.parts_sizes[from].get().fetch_sub(1, std::memory_order_relaxed);

        td.parts_weights[from] -= this_nodes_weight;
        ++td.parts_sizes[to];
        --td.parts_sizes[from];

        //update gain of neighbors / the boundaries have allready been updated
        forall_out_edges(td.G, e, node) {
                ++td.scaned_neighbours;
                NodeID target = td.G.getEdgeTarget(e);
                PartitionID targets_to;
                EdgeWeight ext_degree; // the local external degree

                if (queue->contains(target)) {
                        PartitionID target_from = td.get_local_partition(target);

                        Gain gain = td.compute_gain(target, target_from, targets_to, ext_degree);

                        if (td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
                                return false;
                        }

                        assert(td.moved_idx[target].load(std::memory_order_relaxed));
                        if (ext_degree > 0) {
                                queue->changeKey(target, gain);
                                td.set_local_partition_to_move(target, targets_to);
                        } else {
                                queue->deleteNode(target);
                                td.remove_local_partition_to_move(target);
                        }
                } else {
                        // target was removed from priority queue
                        if (td.moved_idx[target].load(std::memory_order_relaxed)) {
                                continue;
                        }

                        PartitionID target_from = td.get_local_partition(target);

                        Gain gain = td.compute_gain(target, target_from, targets_to, ext_degree);

                        if (td.num_threads_finished.load(std::memory_order_acq_rel) > 0) {
                                return false;
                        }

                        if (ext_degree > 0) {
                                bool expected = false;
                                if (td.moved_idx[target].compare_exchange_strong(expected, true,
                                                                                 std::memory_order_relaxed)) {
                                        queue->insert(target, gain);
                                        td.set_local_partition_to_move(target, targets_to);
                                        td.moved.push_back(target);
                                }
                        }
                }
        } endfor

        return true;
}
}
