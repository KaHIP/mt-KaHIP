#include "data_structure/parallel/algorithm.h"
#include "data_structure/graph_access.h"
#include "data_structure/parallel/thread_pool.h"
#include "data_structure/parallel/time.h"

#include <tbb/concurrent_queue.h>

#include <vector>

namespace parallel {

class numa_aware_graph;

class numa_aware_graph_handle {
public:
        inline NodeID getEdgeTarget(EdgeID edge) {
                return m_edges[edge];
        }

        inline PartitionID getPartitionIndex(NodeID node) {
                return m_G.getPartitionIndex(node);
        }

        inline EdgeID number_of_edges() {
                return m_G.number_of_edges();
        }

        inline NodeID number_of_nodes() {
                return m_G.number_of_nodes();
        }

        inline EdgeID get_first_edge(NodeID node) {
                return m_G.get_first_edge(node);
        }

        inline EdgeID get_first_invalid_edge(NodeID node) {
                return m_G.get_first_invalid_edge(node);
        }

        inline NodeWeight getNodeWeight(NodeID node) {
                return m_G.getNodeWeight(node);
        }

        inline EdgeWeight getNodeDegree(NodeID node) {
                return m_G.getNodeDegree(node);
        }

        friend numa_aware_graph;
private:
        numa_aware_graph_handle(graph_access& G, EdgeID* edges)
                :       m_G(G)
                ,       m_edges(edges)
        {}

        graph_access& m_G;
        EdgeID* m_edges;
};

class numa_aware_graph {
public:
        using handle_type = numa_aware_graph_handle;

        numa_aware_graph(graph_access& G, uint32_t num_threads, uint32_t thread_per_socket)
                :       m_G(G)
                ,       m_num_threads(num_threads)
                ,       m_threads_per_socket(thread_per_socket)
                ,       m_num_sockets(num_threads / thread_per_socket + (num_threads % thread_per_socket != 0 ? 1 : 0))
                ,       m_edges(m_num_sockets, nullptr)
        {}

        void construct() {
                CLOCK_START;
                parallel::submit_for_all([&, this] (uint32_t thread_id) {
                        if (thread_id % m_threads_per_socket == 0) {
                                m_edges[get_socket_id(thread_id)] = reinterpret_cast<EdgeID*>(::operator new(sizeof(EdgeID) * m_G.number_of_edges()));
                        }
                });
                CLOCK_END("Allocate mem");

                CLOCK_START_N;
                Cvector<std::atomic<uint32_t>> offsets_edges(m_num_sockets);
                parallel::submit_for_all([&, this] (uint32_t thread_id) {
                        auto* edges = m_edges[get_socket_id(thread_id)];
                        uint32_t block_size = (uint32_t) sqrt(m_G.number_of_edges());
                        block_size = std::max(block_size, 1000u);
                        auto& offset_edge = offsets_edges[get_socket_id(thread_id)].get();

                        while (true) {
                                uint32_t begin = offset_edge.fetch_add(block_size, std::memory_order_relaxed);
                                uint32_t end = begin + block_size;
                                end = end <= m_G.number_of_edges() ? end : m_G.number_of_edges();

                                if (begin >= m_G.number_of_edges()) {
                                        break;
                                }

                                for (NodeID e = begin; e != end; ++e) {
                                        edges[e] = m_G.graphref->m_edges[e].target;
                                }
                        }
                });
                handles.reserve(m_num_sockets);
                for (size_t socket_id = 0; socket_id < m_num_sockets; ++socket_id) {
                        handles.push_back(numa_aware_graph_handle(m_G, m_edges[socket_id]));
                }
                CLOCK_END("Copy graph");
        }

        inline handle_type& get_handle(uint32_t thread_id) {
                uint32_t socket_id = get_socket_id(thread_id);
                return handles[socket_id];
        }

        ~numa_aware_graph() {
                CLOCK_START;
                parallel::submit_for_all([&, this] (uint32_t thread_id) {
                        if (thread_id % m_threads_per_socket == 0) {
                                ::operator delete(m_edges[get_socket_id(thread_id)]);
                        }
                });
                CLOCK_END("Delete mem");
        }
private:
        graph_access& m_G;
        const uint32_t m_num_threads;
        const uint32_t m_threads_per_socket;
        const uint32_t m_num_sockets;
        std::vector<EdgeID*> m_edges;
        std::vector<handle_type> handles;

        inline uint32_t get_socket_id(uint32_t thread_id) const {
                return thread_id / m_threads_per_socket;
        }
};

template <typename TGraph, typename TFunctor>
void parallel_for_each_node(TGraph& G, TFunctor&& functor) {
        using block_type = std::vector<NodeID>;
        tbb::concurrent_queue <block_type> queue;

        const size_t edge_block_size = std::max((size_t) sqrt(G.number_of_edges()), 1000ul);

        CLOCK_START;
        std::atomic<size_t> offset(0);
        size_t node_block_size = (size_t) sqrt(G.number_of_nodes());
        node_block_size = std::max(node_block_size, 1000ul);
        parallel::submit_for_all([&](uint32_t) {
                block_type block;
                block.reserve(100);
                size_t cur_block_size = 0;
                while (true) {
                        size_t begin = offset.fetch_add(node_block_size, std::memory_order_relaxed);
                        size_t end = begin + node_block_size;
                        end = end <= G.number_of_nodes() ? end : G.number_of_nodes();

                        if (begin >= G.number_of_nodes()) {
                                break;
                        }

                        for (NodeID node = begin; node != end; ++node) {
                                block.push_back(node);
                                auto node_degree = G.getNodeDegree(node);
                                cur_block_size += node_degree > 0 ? node_degree : 1;
                                if (cur_block_size >= edge_block_size) {
                                        // block is full
                                        queue.push(std::move(block));
                                        block.clear();
                                        block.reserve(100);
                                        cur_block_size = 0;
                                }
                        }
                }
                if (!block.empty()) {
                        queue.push(std::move(block));
                }
        });
        CLOCK_END("Init queue");

        CLOCK_START_N;
        parallel::submit_for_all([&](uint32_t id) {
                block_type cur_block;

                while (queue.try_pop(cur_block)) {
                        for (NodeID node : cur_block) {
                                functor(G, node);
                        }
                }
        });
        CLOCK_END("Process graph");
}

}
