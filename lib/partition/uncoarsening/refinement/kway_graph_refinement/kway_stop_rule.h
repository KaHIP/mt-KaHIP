/******************************************************************************
 * kway_stop_rule.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef KWAY_STOP_RULE_ULPK0ZTF
#define KWAY_STOP_RULE_ULPK0ZTF

#include "partition/partition_config.h"
#include <data_structure/parallel/hash_table.h>

#include <cmath>
#include <ios>
#include <fstream>
#include <deque>
#include <iterator>

#undef OUTPUT_GLOBAL_STAT

class kway_stop_rule {
public:
        kway_stop_rule(PartitionConfig & config) {};
        kway_stop_rule() {};
        virtual ~kway_stop_rule() {};

        virtual void push_statistics(Gain gain) = 0;
        virtual void reset_statistics() = 0;
        virtual bool search_should_stop(unsigned int min_cut_idx, 
                                        unsigned int cur_idx, 
                                        unsigned int search_limit) = 0;

#ifdef OUTPUT_GLOBAL_STAT
        virtual void add_stat(uint32_t movement, Gain improvement) {}
#endif
};

class kway_simple_stop_rule : public kway_stop_rule {
public:
        kway_simple_stop_rule(PartitionConfig & config) {};
        virtual ~kway_simple_stop_rule() {};

        void push_statistics(Gain gain) {};
        void reset_statistics() {};
        bool search_should_stop(unsigned int min_cut_idx, 
                                unsigned int cur_idx, 
                                unsigned int search_limit);
};

inline bool kway_simple_stop_rule::search_should_stop(unsigned int min_cut_idx, 
                                                      unsigned int cur_idx, 
                                                      unsigned int search_limit) {
        return cur_idx - min_cut_idx > search_limit;
}

#define TEST_STOPPING_RULE
#ifndef TEST_STOPPING_RULE
class kway_adaptive_stop_rule : public kway_stop_rule {
public:
        kway_adaptive_stop_rule(PartitionConfig & config) : m_steps(0),
                                                            m_expected_gain(0.0),
                                                            m_expected_variance2(0.0),
                                                            pconfig(&config) {}
        virtual ~kway_adaptive_stop_rule() {};

        void push_statistics(Gain gain) {
                //erwartungstreue schätzer für varianz und erwartungswert
                m_expected_gain *= m_steps;
                m_expected_gain += gain;
                if(m_steps == 0) {
                        m_expected_variance2 = 0.0;
                } else {
                        m_expected_variance2 *= (m_steps-1);
                        //expected_variance2 += (gain - expected_gain)*(gain - expected_gain);
                        //real implementation in kaspar
                        m_expected_variance2 += (gain )*(gain );
                }
                m_steps++;

                m_expected_gain /= m_steps;
                if(m_steps > 1)
                        m_expected_variance2 /= (m_steps-1);
        };

        void reset_statistics() {
                m_steps              = 0;
                m_expected_gain      = 0.0;
                m_expected_variance2 = 0.0;
        };

        bool search_should_stop(unsigned int min_cut_idx,
                                unsigned int cur_idx,
                                unsigned int search_limit);
private:
        unsigned m_steps;
        double   m_expected_gain;
        double   m_expected_variance2;
        PartitionConfig * pconfig;
};

inline bool kway_adaptive_stop_rule::search_should_stop(unsigned int min_cut_idx,
                                                        unsigned int cur_idx,
                                                        unsigned int search_limit) {

        return m_steps*m_expected_gain*m_expected_gain >
               pconfig->kway_adaptive_limits_alpha * m_expected_variance2 + pconfig->kway_adaptive_limits_beta && (m_steps != 1);
}
#else
#undef OUT_ADAPTIVE
class kway_adaptive_stop_rule : public kway_stop_rule {
public:
        kway_adaptive_stop_rule(PartitionConfig & config)
                :       m_steps(0)
                ,       m_expected_gain(0.0)
                ,       m_sum_squared_gain(0.0)
                ,       m_expected_variance2(0.0)
                ,       pconfig(&config)
#ifdef OUT_ADAPTIVE
                ,       ftxt("out_a.log", std::ios_base::app | std::ios_base::out)
#endif
                {
#ifdef OUT_ADAPTIVE
                ftxt << "START" << std::endl;
#endif
                }

        virtual ~kway_adaptive_stop_rule() {};

        void push_statistics(Gain gain) {
                //erwartungstreue schätzer für varianz und erwartungswert
                m_expected_gain *= m_steps;
                m_expected_gain += gain;

                m_sum_squared_gain += gain * gain;

                m_steps++;

                m_expected_gain /= m_steps;
                calc_variance();
#ifdef OUT_ADAPTIVE
                ftxt << "{ 'gain' : " << gain;
#endif
        };

        void reset_statistics() {
                m_steps              = 0;
                m_expected_gain      = 0.0;
                m_sum_squared_gain   = 0.0;
                m_expected_variance2 = 0.0;
        };

        bool search_should_stop(unsigned int min_cut_idx,
                                unsigned int cur_idx,
                                unsigned int search_limit);

#ifdef OUTPUT_GLOBAL_STAT
        static void dump_statistics() {
                ALWAYS_ASSERT(m_stat_movements.size() == m_stat_gains.size());
                std::ofstream ftxt("out_a.log");
                if (m_stat_movements.empty()) {
                        return;
                }

                for (size_t i = 0; i < m_stat_movements.size(); ++i) {
                        ftxt << m_stat_movements[i] << '\n';
                }
                ftxt << "end" << std::endl;

                for (size_t i = 0; i < m_stat_gains.size(); ++i) {
                        ftxt << m_stat_gains[i] << '\n';
                }
                ftxt << "end" << std::endl;
        }

        virtual void add_stat(uint32_t movement, Gain improvement) {
                m_stat_movements.push_back(movement);
                m_stat_gains.push_back(improvement);
        }

public:
        static std::vector<uint32_t> m_stat_movements;
        static std::vector<Gain> m_stat_gains;

#endif
private:
        unsigned m_steps;
        double   m_expected_gain;
        Gain   m_sum_squared_gain;
        double   m_expected_variance2;
        PartitionConfig * pconfig;
#ifdef OUT_ADAPTIVE
        std::ofstream ftxt;
#endif

        void calc_variance() {
                if (m_steps > 1) {
                        Gain m_sum_gain = m_expected_gain * m_steps;
                        m_expected_variance2 = (m_sum_squared_gain - 2 * m_expected_gain * m_sum_gain +
                                m_expected_gain * m_expected_gain * m_steps) / (m_steps - 1.0);
                } else {
                        m_expected_variance2 = 0.0;
                }
        }
};

inline bool kway_adaptive_stop_rule::search_should_stop(unsigned int min_cut_idx,
                                                        unsigned int cur_idx,
                                                        unsigned int search_limit) {
#ifdef OUT_ADAPTIVE
        if (m_steps != 1) {
                ftxt << ", 'E' : " << m_expected_gain
                     << ", 'Var' : " << m_expected_variance2
                     << ", 'p * E * E' : " << m_steps*m_expected_gain*m_expected_gain
                     << ", 'a * Var + b' : " << pconfig->kway_adaptive_limits_alpha * m_expected_variance2 + pconfig->kway_adaptive_limits_beta
                     << std::endl;
        }
#endif
        return m_steps*m_expected_gain*m_expected_gain >
                pconfig->kway_adaptive_limits_alpha * m_expected_variance2 + pconfig->kway_adaptive_limits_beta && (m_steps != 1);
}
#endif

class kway_chebyshev_adaptive_stop_rule : public kway_stop_rule {
public:
        kway_chebyshev_adaptive_stop_rule(PartitionConfig& config)
                :       m_steps(0)
                ,       m_total_gain(0)
                ,       m_max_gain(1)
                ,       m_gain_expectation_estimation(0.0)
                ,       probability(0.0)
                ,       m_config(config)
                //,       ftxt("out.log", std::ios_base::app | std::ios_base::out)
        {
                //ftxt << "START" << std::endl;
        }

        virtual ~kway_chebyshev_adaptive_stop_rule() {
        }

        void push_statistics(Gain gain) {
                m_max_gain = std::max(m_max_gain, gain);

                m_total_gain += gain;
                //m_total_gain_squares += gain * gain;

                m_gain_expectation_estimation *= m_steps;
                probability *= m_steps;

                m_gain_expectation_estimation += gain;

                const double t = 1.0 / m_max_gain;
                // E[exp(t*gain)]/E[exp(t)]
                probability += (std::exp(t * (gain - 1)));

                ++m_steps;
                m_gain_expectation_estimation /= m_steps;
                probability /= m_steps;
        }

        void reset_statistics() {
                m_total_gain = 0;
                m_steps = 0;
                m_gain_expectation_estimation = 0.0;
                probability = 0.0;
                //ftxt << "RESET STATISTICS" << std::endl;
        }

        bool search_should_stop(uint32_t min_cut_idx, uint32_t cur_idx, uint32_t search_limit) {
                if (m_steps > 1 && m_gain_expectation_estimation <= 0.0) {
                        const double p = probability;

                        return std::pow(p, int(1 - m_total_gain)) < m_stop_probability;
                } else {
                        return false;
                }
        }
private:
        static constexpr double m_stop_probability = 0.1;
        //std::ofstream // ftxt;

        uint32_t m_steps;
        Gain m_total_gain;
        Gain m_max_gain;
        double m_gain_expectation_estimation;
        double probability;
        const PartitionConfig& m_config;
};

#undef ADAPTIVE_IMPROVED
#define ADAPTIVE

#define START_COMBINED_WITH_ADAPTIVE // best, really helps
#undef START_DEFAULT

#define IMPROVED_DISTRIBUTION

#undef NON_CONST_PROBABILITY
#undef USE_DEQUE

#define OUPUT
class kway_chernoff_adaptive_stop_rule : public kway_stop_rule {
public:
        kway_chernoff_adaptive_stop_rule(PartitionConfig& config)
                :       m_stop_probability(config.chernoff_stop_probability)
                ,       m_gradient_descent_num_steps(config.chernoff_gradient_descent_num_steps)
                ,       m_gradient_descent_step_size(config.chernoff_gradient_descent_step_size)
                ,       m_min_step_limit(config.chernoff_min_step_limit)
                ,       m_max_step_limit(config.chernoff_max_step_limit)
                ,       m_step_limit(0)
                ,       m_steps(0)
                ,       m_total_gain(0)
                ,       m_max_gain(1)
                ,       m_sum_squared_gain(0)
                ,       m_expected_variance2(0.0)
                ,       m_t(1.0)
                ,       m_first(true)
#ifndef USE_DEQUE
                ,       m_gains(32)
#endif
                ,       m_config(config)
                ,       m_adaptive_stop_rule(config)
#ifdef OUPUT
                ,       ftxt("out.log", std::ios_base::app | std::ios_base::out)
#endif
        {
#ifdef OUPUT
                ftxt << "START" << std::endl;
#endif
        }

        virtual ~kway_chernoff_adaptive_stop_rule() {
        }

        static std::string get_algo_name() {
                std::string name = "chernoff_adaptive";
#ifdef ADAPTIVE_IMPROVED
                name += "_improved";
#endif
#ifdef USE_DEQUE
                name += "_deque";
#endif
#ifdef IMPROVED_DISTRIBUTION
                name += "_improved_distr1";
#endif
                return name;
        }
        
        void push_statistics(Gain gain) {
#ifdef OUPUT
                ftxt << "{'gain' : " << gain
                     << ", 'steps': " << m_steps;
#endif

                ++m_steps;
                m_total_gain += gain;
                m_max_gain = std::max(m_max_gain, gain);
#ifdef USE_DEQUE
                m_gains.push_back(gain);
                if (m_gains.size() == m_max_step_limit) {
                        m_gains.pop_front();
                }
#else
                ++m_gains[gain];
#endif

#ifdef START_COMBINED_WITH_ADAPTIVE
                m_adaptive_stop_rule.push_statistics(gain);
#endif
#ifdef IMPROVED_DISTRIBUTION
                m_sum_squared_gain += gain * gain;
#endif
        }

        void reset_statistics() {
                reset_adaptive_statistics();
        }

        bool search_should_stop(uint32_t min_cut_idx, uint32_t cur_idx, uint32_t search_limit) {
#ifdef START_DEFAULT
                bool res = false;
                if (m_steps >= m_min_step_limit && m_total_gain <= 0.0) {
                        if (m_first) {
                                m_first = false;
                                m_t = 1.0 / m_max_gain;
                        }
#ifdef IMPROVED_DISTRIBUTION
                        calc_variance();
#endif
                        const double p = get_probability();

                        double cur_stop_probability = m_stop_probability;

                        res = p < cur_stop_probability;

#ifdef OUPUT
                        ftxt << ", t : " << m_t
                             << ", n : " << get_n()
                             << ", p : " << p
                             << ", stop_p: " << cur_stop_probability
                             << ", total_gain: " << m_total_gain
                             << "}," << std::endl;
                        ftxt << "m_gains = {";
                        for (const auto& data : m_gains) {
                                ftxt << data.first << " : " << (data.second + 0.0) / m_steps << ", ";
                        }
                        ftxt << "}" << std::endl;
#endif
                }
#endif


#ifdef START_COMBINED_WITH_ADAPTIVE
                bool res = false;
                if (m_steps >= m_min_step_limit) {
                        if (m_total_gain <= 0.0) {
                                if (m_first) {
                                        m_first = false;
                                        m_t = 1.0 / m_max_gain;
                                }
#ifdef IMPROVED_DISTRIBUTION
                                calc_variance();
#endif
                                const double p = get_probability();

                                double cur_stop_probability = m_stop_probability;

                                res = p < cur_stop_probability;

#ifdef OUPUT
                                double expectation = get_expectation();
                                double std_dev = sqrt(m_expected_variance2);

                                ftxt << ", 't' : " << m_t
                                     << ", 'p' : " << p
                                     << ", 'stop_p': " << cur_stop_probability
                                     << ", 'total_gain': " << m_total_gain
                                     << ", 'step': " << m_steps
                                     << ", 'E' : " << expectation
                                     << ", 'std': " << std_dev
                                     << "}," << std::endl;
                                ftxt << "m_gains = {";
                                for (const auto& data : m_gains) {
                                        ftxt << data.first << " : " << (data.second + 0.0) / m_steps << ", ";
                                }
                                ftxt << "}" << std::endl;
#endif
                        }
                } else {
                        res = m_adaptive_stop_rule.search_should_stop(min_cut_idx, cur_idx, search_limit);
                }
#endif

#ifndef USE_DEQUE
                if (m_steps == m_max_step_limit) {
                        reset_chernoff_statistics();
                }
#endif
                return res;
        }

#ifdef OUTPUT_GLOBAL_STAT
        static void dump_statistics() {
                ALWAYS_ASSERT(m_stat_movements.size() == m_stat_gains.size());
                std::ofstream ftxt("out_a.log");
                if (m_stat_movements.empty()) {
                        return;
                }

                for (size_t i = 0; i < m_stat_movements.size(); ++i) {
                        ftxt << m_stat_movements[i] << '\n';
                }
                ftxt << "end" << std::endl;

                for (size_t i = 0; i < m_stat_gains.size(); ++i) {
                        ftxt << m_stat_gains[i] << '\n';
                }
                ftxt << "end" << std::endl;
        }

        virtual void add_stat(uint32_t movement, Gain improvement) {
                m_stat_movements.push_back(movement);
                m_stat_gains.push_back(improvement);
        }

public:
        static std::vector<uint32_t> m_stat_movements;
        static std::vector<Gain> m_stat_gains;

#endif
private:
        double m_stop_probability;
        uint32_t m_gradient_descent_num_steps;
        double m_gradient_descent_step_size;
        uint32_t m_min_step_limit;
        uint32_t m_max_step_limit;
        uint32_t m_step_limit;

        uint32_t m_steps;
        Gain m_total_gain;
        Gain m_max_gain;
#ifdef IMPROVED_DISTRIBUTION
        Gain m_sum_squared_gain;
        double m_expected_variance2;
#endif
        double m_t;
        bool m_first;
#ifdef USE_DEQUE
        std::deque<Gain> m_gains;
#else
        parallel::hash_map<Gain, uint32_t> m_gains;
#endif

        const PartitionConfig& m_config;
        kway_adaptive_stop_rule m_adaptive_stop_rule;

#ifdef OUPUT
        mutable std::ofstream ftxt;
#endif

        void reset_adaptive_statistics() {
#ifdef START_COMBINED_WITH_ADAPTIVE
                m_adaptive_stop_rule.reset_statistics();
#endif
        }

        void reset_chernoff_statistics() {
                m_steps = 0;
                m_total_gain = 0;
                m_max_gain = 1;
                m_t = 1.0;
                m_first = true;
                m_gains.clear();
                m_step_limit = 0;
                m_sum_squared_gain = 0;
                m_expected_variance2 = 0.0;
        }

        template <typename function_type>
        double get_expectation(function_type function) const {
                double expectation = 0;
                for (const auto& data : m_gains) {
#ifdef USE_DEQUE
                        expectation += (function(data) + 0.0) / m_gains.size();
#else
                        expectation += (data.second + 0.0) / m_steps * function(data.first);
#endif
                }
                return expectation;
        }

        template <typename function_type, typename predicate_type>
        double get_partial_expectation(function_type function, predicate_type condition) const {
                double expectation = 0;
                for (const auto& data : m_gains) {
                        if (condition(data.first)) {
                                expectation += (data.second + 0.0) / m_steps * function(data.first);
                        }
                }
                return expectation;
        }

        template <typename function_type, typename predicate_type>
        double get_conditional_expectation(function_type function, predicate_type condition) const {
                double expectation = 0;

                uint32_t num_condition_true = 0;
                for (const auto& data : m_gains) {
#ifdef USE_DEQUE
                        if (condition(data)) {
                                ++num_condition_true;
                        }
#else
                        if (condition(data.first)) {
                                num_condition_true += data.second;
                        }
#endif
                }

                for (const auto& data : m_gains) {
#ifdef USE_DEQUE
                        if (condition(data)) {
                                expectation += (function(data) + 0.0) / num_condition_true;
                        }
#else
                        if (condition(data.first)) {
                                expectation += (data.second + 0.0) / num_condition_true * function(data.first);
                        }
#endif
                }
                return expectation;
        }

        // X = sum_i X_i where X_i is random variable -- gain. Using chernoff bound
        // we receive that P[X >= 1] <= min_t E[exp(tX)]/exp(t). We assume that all X_i are independent and
        // then E[exp(tX)] = E[exp(tX_1)*exp(tX_2)*...*exp(tX_n)] = E[exp(t*X_1)] * E[exp(t*X_2)] * ... * E[exp(t*X_n)] =
        //  E[exp(t*X_1)] ^ n
        double get_probability() {
#ifdef ADAPTIVE_IMPROVED
                int n = 0;
                if (m_step_limit == 0) {
                        n = get_n();

                        if (n == -1) {
                                // there is no any positive gain
                                return 0.0;
                        }
                        m_step_limit = n;
                }
                n = m_step_limit;
                --m_step_limit;
#endif
#ifdef ADAPTIVE
                int n = get_n();

                if (n == -1) {
                        // there is no any positive gain
                        return 0.0;
                }

                n = n > 0 ? n : 10;

#endif

#ifdef OUPUT
                ftxt << ", 'n' : " << n;
#endif
                int a = 1 - m_total_gain;
                double cur_val = get_probability(m_t, n, a);
                double min_val = cur_val;
                double min_t = m_t;
                for (size_t i = 0; i < m_gradient_descent_num_steps; ++i) {
                        double gradient_step = m_gradient_descent_step_size;
                        double derivative = get_derivative(m_t, n, a);
                        double new_t = m_t - gradient_step * derivative;

                        double new_val = get_probability(new_t, n, a);
                        while (i < m_gradient_descent_num_steps && (new_t <= 0.0 || new_val > cur_val)) {
                                gradient_step *= 0.1;
                                new_t = m_t - gradient_step * derivative;
                                new_val = get_probability(new_t, n, a);

                                ++i;
                        }

                        if (new_t > 0.0 && new_val < cur_val) {
                                m_t = new_t;
                                cur_val = new_val;
                        }

                        if (cur_val < min_val) {
                                min_val = cur_val;
                                min_t = m_t;
                        }
                }

                m_t = min_t;

                // E[exp(t * gain)] ^ n / exp(t)
                return min_val;
        }

        double get_probability(double t, int n, int a) const {
#ifdef IMPROVED_DISTRIBUTION
                double expectation = get_expectation();
                double std_dev = sqrt(m_expected_variance2);

                //n = expectation != 0.0 ? m_expected_variance2 / (4.0 * expectation * expectation) : 10;
                double f = get_conditional_expectation([t](Gain gain) {
                        return std::exp(t * gain);
                }, [=] (Gain gain) {
                        return exp_cond(gain, expectation, std_dev);
                }
                );
#else
                double f = get_expectation([t, a](Gain gain) {
                        return std::exp(t * gain);
                });
#endif
                return std::pow(f, n) / std::exp(t * a);
        }

#ifdef IMPROVED_DISTRIBUTION

        inline bool exp_cond(Gain gain, double expectation, double std_dev) const {
                //return gain >= expectation && gain <= expectation + 3 * std_dev;
                return gain >= expectation;
        }

        inline double get_expectation() const {
                return m_total_gain / (m_steps + 0.0);
        }

        void calc_variance() {
                if (m_steps > 1) {
                        double expected_gain = (m_total_gain + 0.0) / m_steps;
                        m_expected_variance2 = (m_sum_squared_gain - 2 * expected_gain * m_total_gain +
                                expected_gain * expected_gain * m_steps) / (m_steps - 1.0);
                } else {
                        m_expected_variance2 = 0.0;
                }
        }
#endif

        double get_derivative(double t, int n, int a) const {
#ifdef IMPROVED_DISTRIBUTION
                double expectation = get_expectation();
                double std_dev = sqrt(m_expected_variance2);

                double f1 = get_conditional_expectation([t](Gain gain) {
                        return std::exp(t * gain);
                }, [=] (Gain gain) {
                        return exp_cond(gain, expectation, std_dev);
                });

                double f2 = get_conditional_expectation([t, n, a](Gain gain) {
                        return std::exp(t * gain) * (n * gain - a);
                }, [=] (Gain gain) {
                        return exp_cond(gain, expectation, std_dev);
                });
#else
                double f1 = get_expectation([t](Gain gain) {
                        return std::exp(t * gain);
                });

                double f2 = get_expectation([t, n, a](Gain gain) {
                        return std::exp(t * gain) * (n * gain - a);
                });
#endif
                return std::pow(f1, n - 1) * f2 / std::exp(t * a);
        }

        int get_n() const {
#ifdef IMPROVED_DISTRIBUTION
                // m_total_gain <= 0
                // int(1 - m_total_gain) / E[X_1|X_1 >= 1]
//                double expectation = get_expectation();
//                double std_dev = sqrt(m_expected_variance2);

//                double cond_expect = get_conditional_expectation([](Gain gain) {
//                        return gain;
//                }, [=](Gain gain) {
//                        //return exp_cond(gain, expectation, std_dev);
//                        //return gain >= expectation && gain <= expectation + std_dev;
//                        return gain >= 1;
//                });
                double cond_expect = get_partial_expectation([](Gain gain) {
                        return gain;
                }, [=](Gain gain) {
                        //return exp_cond(gain, expectation, std_dev);
                        //return gain >= expectation && gain <= expectation + std_dev;
                        return gain >= 1;
                });
#else
                // m_total_gain <= 0
                // int(1 - m_total_gain) / E[X_1|X_1 >= 1]
                double cond_expect = get_conditional_expectation([](Gain gain) {
                        return gain;
                }, [](Gain gain) {
                        return gain >= 1;
                });
#endif
                return cond_expect > 0.0 ? (1 - m_total_gain) / cond_expect : -1;
        }
};

#endif /* end of include guard: KWAY_STOP_RULE_ULPK0ZTF */
