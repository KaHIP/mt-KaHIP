#pragma once

#include "data_structure/parallel/cache.h"
#include "data_structure/parallel/metaprogramming_utils.h"
#include "data_structure/parallel/thread_pool.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <future>
#include <vector>

#include <ips4o/ips4o.hpp>
#include <omp.h>
#include <parallel/algorithm>
#include <parallel/numeric>

namespace parallel {
template<typename _value_type>
class integer_iterator : public std::iterator<std::random_access_iterator_tag, _value_type> {
private:
        using base = std::iterator<std::random_access_iterator_tag, _value_type>;
public:
        typedef typename base::iterator_category iterator_category;
        typedef typename base::value_type value_type;
        typedef typename base::difference_type difference_type;
        typedef typename base::pointer pointer;
        typedef typename base::reference reference;

        typedef integer_iterator<value_type> self;

        inline integer_iterator(value_type value)
                : m_value(value) {}

        inline value_type operator*() const {
                return m_value;
        }

        inline difference_type operator-(self& other) const {
                return m_value - other.m_value;
        }

        inline self operator+(difference_type val) const {
                return self(m_value + val);
        }

        inline self& operator+=(difference_type val) {
                m_value += val;
                return *this;
        }

        inline self& operator++() {
                ++m_value;
                return *this;
        }

        inline self operator++(int) {
                value_type old = m_value;
                ++m_value;
                return self(old);
        }

        inline bool operator==(const self& other) const {
                return m_value == other.m_value;
        }

        inline bool operator!=(const self& other) const {
                return m_value != other.m_value;
        }

private:
        static_assert(std::is_integral<value_type>::value, "Value should be of integral type");

        value_type m_value;
};

template<typename Iterator, typename Function>
Function apply_to_range(Iterator first, Iterator last, Function f)
{
        f(first, last);
        return f;
}

template<bool Sync, template <typename, typename> typename AggregateFuncWrapper,
        typename Iterator, typename Functor, typename TaskConsumer>
static Functor for_each(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func,
                        std::random_access_iterator_tag) {
        size_t size = std::distance(first, last);
        size_t proc = task_consumer.NumThreads() + 1;
        AggregateFuncWrapper<Iterator, Functor> wrapper;

        if (size <= 5 * proc) {
                return wrapper(first, last, func);
        }

        size_t work_per_thread = size / proc;

        std::vector<std::future<Functor>> futures;
        if (Sync) {
                futures.reserve(proc - 1);
        }

        for (size_t i = 0; i + 1 < proc; ++i) {
                if (Sync) {
                        futures.push_back(task_consumer.Submit(i, wrapper, first,
                                                               first + work_per_thread, func));
//                        futures.push_back(task_consumer.Submit(wrapper, first,
//                                                               first + work_per_thread, func));
                } else {
                        task_consumer.Submit(i, wrapper,
                                             first, first + work_per_thread, func);
//                        task_consumer.Submit(wrapper,
//                                             first, first + work_per_thread, func);
                }
                first += work_per_thread;
        }

        wrapper(first, last, func);

        if (Sync) {
                std::for_each(futures.begin(), futures.end(), [](auto& future) {
                        future.get();
                });
        }
        return func;
}

template <typename iterator_type, bool is_integral = std::is_integral<iterator_type>::value>
struct iterator
{
        using type = iterator_type;
        static type& get(type& iter) {
                return iter;
        }
};

template <typename value_type>
struct iterator<value_type, true>
{
        using type = integer_iterator<value_type>;

        static type get(value_type val) {
                return integer_iterator<value_type>(val);
        }
};

template <typename Iterator, typename Functor>
struct for_each_wrapper {
        inline Functor operator() (Iterator begin, Iterator end, Functor f) const {
                return std::for_each(begin, end, f);
        }
};

template <typename Iterator, typename Functor>
struct apply_to_range_wrapper {
        inline Functor operator() (Iterator begin, Iterator end, Functor f) const {
                return parallel::apply_to_range(begin, end, f);
        }
};

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor for_each_sync(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;

        return for_each<true, for_each_wrapper>(getter(first), getter(last), task_consumer, func,
                                                typename std::iterator_traits<iter_type>::iterator_category());
}

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor for_each_async(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;
        return for_each<false, for_each_wrapper>(getter(first), getter(last), task_consumer, func,
                                                 typename std::iterator_traits<iter_type>::iterator_category());
}

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor apply_to_range_sync(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;

        return for_each<true, apply_to_range_wrapper>(getter(first), getter(last), task_consumer, func,
                                                      typename std::iterator_traits<iter_type>::iterator_category());
}

template<typename Iterator, typename Functor, typename TaskConsumer>
Functor apply_to_range_async(Iterator first, Iterator last, TaskConsumer& task_consumer, Functor func) {
        using iter_type = typename iterator<Iterator>::type;
        auto getter = iterator<Iterator>::get;
        return for_each<false, apply_to_range_wrapper>(getter(first), getter(last), task_consumer, func,
                                                       typename std::iterator_traits<iter_type>::iterator_category());
}

template <typename T>
using Cvector = std::vector<CacheAlignedData<T>>;

template <typename T>
T* allocate_raw_memory(size_t size) {
        return reinterpret_cast<T*>(sizeof(T) * size);
}

template <typename T>
class ParallelVector {
public:
        using Self = ParallelVector<T>;
        using Type = T;

        explicit ParallelVector(size_t size)
                :       m_ptr(nullptr)
                ,       m_size(size)
        {
                m_ptr = reinterpret_cast<T*>(::operator new(size * sizeof(T)));
        }

        ~ParallelVector() {
                ::operator delete(m_ptr);
        }

        const Type& operator[] (size_t index) const {
                return m_ptr[index];
        }

        Type& operator[] (size_t index) {
                return m_ptr[index];
        }

        const Type& front() const {
                return m_ptr[0];
        }

        Type& front() {
                return m_ptr[0];
        }

        const Type& back() const {
                return m_ptr[size() - 1];
        }

        Type& back() {
                return m_ptr[size() - 1];
        }

        size_t size() const {
                return m_size;
        }

        void swap(Self& other) {
                std::swap(m_ptr, other.m_ptr);
                std::swap(m_size, other.m_size);
        }

        const Type* begin() const {
                return m_ptr;
        }

        Type* begin() {
                return m_ptr;
        }

        const Type* end() const {
                return m_ptr + m_size;
        }

        Type* end() {
                return m_ptr + m_size;
        }

private:
        Type* m_ptr;
        size_t m_size;
};

template<typename Iterator, typename Functor>
static void parallel_for_each(Iterator begin, Iterator end, Functor functor) {
        std::vector<std::future<void>> futures;
        futures.reserve(g_thread_pool.NumThreads());

        std::atomic<size_t> offset(0);
        size_t size = end - begin;
        size_t block_size = (size_t) sqrt(size);
        block_size = std::max(block_size, size_t(1000));
        auto task = [&] (uint32_t thread_id) {
                size_t cur_begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                while (cur_begin < size) {
                        size_t cur_end = std::min(cur_begin + block_size, size);

                        for (Iterator elem = begin + cur_begin; elem != begin + cur_end; ++elem) {
                                if constexpr (function_traits<Functor>::arity == 1) {
                                        functor(*elem);
                                }
                                if constexpr (function_traits<Functor>::arity == 2) {
                                        functor(*elem, thread_id);
                                }
                        }
                        cur_begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }
        };

        for (uint32_t i = 0; i < g_thread_pool.NumThreads(); ++i) {
                futures.push_back(g_thread_pool.Submit(i, task, i + 1));
        }
        task(uint32_t(0));

        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                future.get();
        });
}

template<typename Integer_type, typename Functor>
static void parallel_for_index(Integer_type begin, Integer_type end, Functor functor) {
        static_assert(std::is_integral<Integer_type>::value, "Integral required.");

        std::vector<std::future<void>> futures;
        futures.reserve(g_thread_pool.NumThreads());

        std::atomic<size_t> offset(0);
        size_t size = end - begin;
        size_t block_size = (size_t) sqrt(size);
        block_size = std::max(block_size, size_t(1000));
        auto task = [&] (uint32_t thread_id) {
                size_t cur_begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                while (cur_begin < size) {
                        size_t cur_end = std::min(cur_begin + block_size, size);

                        for (Integer_type elem = cur_begin; elem != cur_end; ++elem) {
                                if constexpr (function_traits<Functor>::arity == 1) {
                                        functor(elem);
                                }
                                if constexpr (function_traits<Functor>::arity == 2) {
                                        functor(elem, thread_id);
                                }
                        }
                        cur_begin = offset.fetch_add(block_size, std::memory_order_relaxed);
                }
        };

        for (uint32_t i = 0; i < g_thread_pool.NumThreads(); ++i) {
                futures.push_back(g_thread_pool.Submit(i, task, i + 1));
        }
        task(0);

        std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                future.get();
        });
}

template <typename Iterator>
void random_shuffle(Iterator begin, Iterator end, uint32_t num_threads) {
        ALWAYS_ASSERT(num_threads > 0);

        parallel::g_thread_pool.Clear();
        parallel::Unpin();
        omp_set_dynamic(false);
        omp_set_num_threads(num_threads);

        __gnu_parallel::random_shuffle(begin, end);

        omp_set_num_threads(0);
        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(num_threads - 1);
}

template <typename InputIterator, typename OutputIterator>
void partial_sum(InputIterator begin, InputIterator end, OutputIterator out, uint32_t num_threads) {
        ALWAYS_ASSERT(num_threads > 0);

        parallel::g_thread_pool.Clear();
        parallel::Unpin();
        omp_set_dynamic(false);
        omp_set_num_threads(num_threads);

        __gnu_parallel::partial_sum(begin, end, out);

        omp_set_num_threads(0);
        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(num_threads - 1);
}

// calculates prefix sum for each prefix not including last element of the prefix
template <typename InputIterator, typename OutputIterator>
void partial_sum_open_interval(InputIterator begin, InputIterator end, OutputIterator out, uint32_t num_threads) {
        ALWAYS_ASSERT(num_threads > 0);

        size_t size = end - begin;

        // copy array
        using value_type = typename InputIterator::value_type;
        ParallelVector<value_type> copy(size);

        // calculate prefix sum for copy
        partial_sum(begin, end, copy.begin(), num_threads);

        // change the input array
        parallel::parallel_for_index(size_t(0), size, [&](size_t index) {
                *(out + index) = copy[index] - *(begin + index);
        });
}

template <typename Iterator, typename Functor>
void sort(Iterator begin, Iterator end, Functor functor, uint32_t num_threads) {
        ALWAYS_ASSERT(num_threads > 0);

        parallel::g_thread_pool.Clear();
        parallel::Unpin();

        ips4o::parallel::sort(begin, end, functor, num_threads);

        parallel::PinToCore(0);
        parallel::g_thread_pool.Resize(num_threads - 1);
}

}