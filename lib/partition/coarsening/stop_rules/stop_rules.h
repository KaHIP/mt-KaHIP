/******************************************************************************
 * stop_rules.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef STOP_RULES_SZ45JQS6
#define STOP_RULES_SZ45JQS6

#include <cmath>
#include <functional>
#include <memory>

#include "data_structure/parallel/cache.h"
#include "partition_config.h"

class stop_rule {
        public:
                stop_rule() {};
                virtual ~stop_rule() {};
                virtual bool stop( NodeID number_of_finer_vertices, graph_access& coarser) = 0;
};

class separator_simple_stop_rule : public stop_rule {
        public:
                separator_simple_stop_rule(PartitionConfig & config, NodeID number_of_nodes) {
	                num_stop = config.sep_num_vert_stop;
                        if(config.disable_max_vertex_weight_constraint) {
                                config.max_vertex_weight = config.upper_bound_partition; 
                        } else {
                                config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                        }
                };

                virtual ~separator_simple_stop_rule() {};
                bool stop( NodeID number_of_finer_vertices, graph_access& coarser);

        private:
                NodeID num_stop;
};

inline bool separator_simple_stop_rule::stop(NodeID no_of_finer_vertices, graph_access& coarser) {
        NodeID no_of_coarser_vertices = coarser.number_of_nodes();
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}

class simple_stop_rule : public stop_rule {
        public:
                simple_stop_rule(PartitionConfig & config, NodeID number_of_nodes) {
                        double x = 60;
	                num_stop = std::max(number_of_nodes/(2.0*x*config.k), 60.0*config.k);
                        if(config.disable_max_vertex_weight_constraint) {
                                config.max_vertex_weight = config.upper_bound_partition; 
                        } else {
                                config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                        }
                };
                virtual ~simple_stop_rule() {};
                bool stop( NodeID number_of_finer_vertices, graph_access& coarser);

        private:
                NodeID num_stop;
};

inline bool simple_stop_rule::stop(NodeID no_of_finer_vertices, graph_access& coarser) {
        NodeID no_of_coarser_vertices = coarser.number_of_nodes();
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}
class strong_stop_rule : public stop_rule {
        public:
                strong_stop_rule(PartitionConfig & config, NodeID number_of_nodes) {
                        num_stop = config.k;
                        config.max_vertex_weight = config.upper_bound_partition;
                };
                virtual ~strong_stop_rule() {};
                bool stop( NodeID number_of_finer_vertices, graph_access& coarser);

        private:
                NodeID num_stop;
};

inline bool strong_stop_rule::stop(NodeID no_of_finer_vertices, graph_access& coarser) {
        NodeID no_of_coarser_vertices = coarser.number_of_nodes();
        double contraction_rate = 1.0 * no_of_finer_vertices / (double)no_of_coarser_vertices;
        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
}

class multiple_k_stop_rule : public stop_rule {
        public:
                explicit multiple_k_stop_rule (PartitionConfig & config, NodeID number_of_nodes = 0, EdgeID number_of_edges = 0) {
                        num_stop = config.num_vert_stop_factor*config.k;

                        if(config.disable_max_vertex_weight_constraint) {
                                config.max_vertex_weight = config.upper_bound_partition; 
                        } else {
                                if(config.initial_partitioning) {
                                        //if we perform initial partitioning we relax this constraint
                                        config.max_vertex_weight = 1.5*((double)config.work_load)/(2*config.num_vert_stop_factor); 
                                } else {
                                        config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                                }
                        }

                };
                virtual ~multiple_k_stop_rule () {};
                bool stop( NodeID number_of_finer_vertices, graph_access& coarser) override {
                        NodeID no_of_coarser_vertices = coarser.number_of_nodes();
                        double contraction_rate = 1.0 * number_of_finer_vertices / (double)no_of_coarser_vertices;
                        return contraction_rate >= 1.1 && no_of_coarser_vertices >= num_stop;
                }

                void reset() {
                }
        private:
                NodeID num_stop;
};

class mem_stop_rule : public stop_rule {
public:
        mem_stop_rule (PartitionConfig & config, NodeID number_of_nodes)
                :       m_config(config)
        {
                NodeID num_stop = config.num_vert_stop_factor*config.k;

                if(config.disable_max_vertex_weight_constraint) {
                        config.max_vertex_weight = config.upper_bound_partition;
                } else {
                        if(config.initial_partitioning) {
                                //if we perform initial partitioning we relax this constraint
                                config.max_vertex_weight = 1.5*((double)config.work_load)/(2*config.num_vert_stop_factor);
                        } else {
                                config.max_vertex_weight = (NodeWeight)(1.5*config.work_load/num_stop);
                        }
                }

        };
        virtual ~mem_stop_rule () {};
        bool stop(NodeID, graph_access& coarser) override {
                return coarser.mem() > m_config.l3_cache_size / 2;
        }

private:
        const PartitionConfig& m_config;
};

class multiple_k_strong_contraction : public multiple_k_stop_rule {
public:
        multiple_k_strong_contraction(PartitionConfig& _config, NodeID _number_of_nodes, EdgeID _number_of_edges)
                :       multiple_k_stop_rule(_config, _number_of_nodes, _number_of_edges)
                ,       config(_config)
                ,       attemps_to_contract_more(0)
                ,       number_of_nodes(_number_of_nodes)
                ,       number_of_edges(_number_of_edges)
        {}

        virtual ~multiple_k_strong_contraction(){
        }

        bool stop(NodeID number_of_finer_vertices, graph_access& coarser) override {
                bool res = multiple_k_stop_rule::stop(number_of_finer_vertices, coarser);

                if (!res && attemps_to_contract_more < max_attempts &&
                    //((coarser.number_of_edges() + 0.0) / coarser.number_of_nodes() > 1.5 || coarser.number_of_edges() > 500000)
                    //(coarser.number_of_edges() > number_of_edges * 0.001 && coarser.number_of_nodes() > 1000 * config.k)
                    //((coarser.number_of_edges() > number_of_edges * 0.001 || coarser.number_of_edges() > 300000) && coarser.number_of_nodes() > 1000 * config.k)
                    (coarser.number_of_edges() > number_of_edges * 0.001 || coarser.number_of_edges() > 300000)
                ) {
                        config.cluster_coarsening_factor /= 2;
                        config.cluster_coarsening_factor = std::max(config.cluster_coarsening_factor, 1);
                        //config.accept_small_coarser_graphs = false;

                        res = true;
                        ++attemps_to_contract_more;
                }

                return res;
        }

        void reset() {
                attemps_to_contract_more = 0;
        }
private:
        static constexpr uint32_t max_attempts = 2;
        PartitionConfig& config;
        uint32_t attemps_to_contract_more;
        const NodeID number_of_nodes;
        const EdgeID number_of_edges;
};

template <typename main_stop_rule_type>
class stop_rule_with_reset : public main_stop_rule_type {
public:
        template <typename functor_type>
        stop_rule_with_reset(PartitionConfig& _config, functor_type&& _update_config, NodeID _number_of_nodes,
                             EdgeID _number_of_edges)
                :       main_stop_rule_type(_config, _number_of_nodes, _number_of_edges)
                ,       config(_config)
                ,       phase(stop_rule_phase::first_phase)
                ,       update_config(std::forward<functor_type>(_update_config))
        {}

        virtual ~stop_rule_with_reset() {
        }

        bool stop(NodeID number_of_finer_vertices, graph_access& coarser) override {
                bool res = main_stop_rule_type::stop(number_of_finer_vertices, coarser);

                if (!res && phase == stop_rule_phase::first_phase) {
                        res = true;
                        phase = stop_rule_phase::second_phase;
                        update_config(config);
                        main_stop_rule_type::reset();
                }

                return res;
        }

private:
        enum class stop_rule_phase {
                first_phase,
                second_phase
        };

        PartitionConfig& config;
        stop_rule_phase phase;
        std::function<void(PartitionConfig&)> update_config;
};

template <typename first_stop_rule_type, typename second_stop_rule>
class double_stop_rule : public stop_rule {
public:
        template <typename functor_type>
        double_stop_rule(PartitionConfig& _config, functor_type&& _update_config, NodeID _number_of_nodes,
                         EdgeID _number_of_edges)
                :       config(_config)
                ,       number_of_nodes(_number_of_nodes)
                ,       number_of_edges(_number_of_edges)
                ,       phase(stop_rule_phase::first_phase)
                ,       update_config(std::forward<functor_type>(_update_config))
                ,       stop_rule_impl(std::make_unique<first_stop_rule_type>(_config, _number_of_nodes, _number_of_edges))
        {}

        virtual ~double_stop_rule() {
        }

        bool stop(NodeID number_of_finer_vertices, graph_access& coarser) override {
                bool res = stop_rule_impl->stop(number_of_finer_vertices, coarser);

                if (!res && phase == stop_rule_phase::first_phase) {
                        res = true;
                        phase = stop_rule_phase::second_phase;
                        update_config(config);
                        stop_rule_impl = std::make_unique<second_stop_rule>(config, number_of_nodes, number_of_edges);
                }

                return res;
        }

private:
        enum class stop_rule_phase {
                first_phase,
                second_phase
        };

        PartitionConfig& config;
        NodeID number_of_nodes;
        EdgeID number_of_edges;
        stop_rule_phase phase;
        std::function<void(PartitionConfig&)> update_config;
        std::unique_ptr<stop_rule> stop_rule_impl;
};

#endif /* end of include guard: STOP_RULES_SZ45JQS6 */
