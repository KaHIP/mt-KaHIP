#pragma once

#include <cmath>
#include <type_traits>

#include <iostream>

namespace parallel {

template <typename T>
class estimator {
public:
        static_assert(std::is_arithmetic<T>::value);

        void put(T value) {
                sum += value;
                sum_squared += value * value;
                ++num_elems;
        }

        double get_expectation() const {
                if (num_elems > 0) {
                        return (sum + 0.0) / num_elems;
                } else {
                        return 0.0;
                }
        }

        double get_variance() const {
                if (num_elems > 1) {
                        double exp = get_expectation();
                        return (sum_squared - 2 * exp * sum + num_elems * exp * exp + 0.0) / (num_elems - 1);
                } else {
                        return 0.0;
                }
        }

        double get_std_deviation() const {
                return sqrt(get_variance());
        }

        size_t get_num_elems() const {
                return num_elems;
        }

        explicit estimator()
                :       sum (T())
                ,       sum_squared(T())
                ,       num_elems(0)
        {}
private:
        T sum;
        T sum_squared;
        size_t num_elems;
};

}
