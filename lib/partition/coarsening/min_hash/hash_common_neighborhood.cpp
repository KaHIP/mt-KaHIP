#include "data_structure/parallel/hash_function.h"
#include "partition/coarsening/min_hash/hash_common_neighborhood.h"

#include <unordered_map>
#include <iostream>

void hash_common_neighborhood::match(const PartitionConfig& config,
                                     graph_access& G,
                                     Matching&,
                                     CoarseMapping& coarse_mapping,
                                     NodeID& no_of_coarse_vertices,
                                     NodePermutationMap&) {
        find_vertices_with_common_neighbors(config, G, coarse_mapping, no_of_coarse_vertices);
}

void hash_common_neighborhood::find_vertices_with_common_neighbors(const PartitionConfig& config,
                                                                   graph_access& G,
                                                                   CoarseMapping& coarse_mapping,
                                                                   NodeID& no_of_coarse_vertices) {
        const NodeID cluster_upperbound = (NodeID) ceil(
                (config.upper_bound_partition + 0.0) / config.cluster_coarsening_factor);

        std::vector<NodeID> cluster_sizes(G.number_of_nodes());
        std::unordered_map<uint64_t, std::vector<NodeID>> buckets;
        parallel::MurmurHash<NodeID> hash(config.seed);

        for (NodeID node = 0; node < G.number_of_nodes(); ++node) {
                uint64_t hash_value = 0;
                //uint64_t hash_value = std::numeric_limits<uint64_t>::max();
                forall_out_edges(G, e, node){
                        NodeID neighbor = G.getEdgeTarget(e);
                        hash_value ^= hash(neighbor);
                        //hash_value = std::min(hash(neighbor), hash_value);
                } endfor


                buckets[hash_value].push_back(node);
                cluster_sizes[coarse_mapping[node]] += G.getNodeWeight(node);
        }

        for (const auto& bucket : buckets) {
                size_t i = 0;
                while (i < bucket.second.size()) {
                        NodeID node = bucket.second[i];
                        NodeID cluster = coarse_mapping[node];

                        while (++i < bucket.second.size()) {
                                NodeID next_node = bucket.second[i];

                                if (G.getNodeWeight(next_node) + cluster_sizes[cluster] <= cluster_upperbound) {
                                        cluster_sizes[coarse_mapping[next_node]] -= G.getNodeWeight(next_node);
                                        coarse_mapping[next_node] = cluster;
                                        cluster_sizes[cluster] += G.getNodeWeight(next_node);
                                } else {
                                        break;
                                }
                        }
                }
        }

        if (coarse_mapping.empty()) {
                return;
        }

        std::vector<NodeWeight> cluster_map(G.number_of_nodes());
        forall_nodes(G, node){
                                PartitionID cur_cluster = coarse_mapping[node];
                                cluster_map[cur_cluster] = 1;
        } endfor

        for (size_t i = 1; i < cluster_map.size(); ++i)
                cluster_map[i] += cluster_map[i - 1];

        forall_nodes(G, node){
                                coarse_mapping[node] = cluster_map[coarse_mapping[node]] - 1;
        } endfor
        no_of_coarse_vertices = cluster_map.back();
}