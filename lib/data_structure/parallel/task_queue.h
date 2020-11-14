#pragma once

#include "data_structure/parallel/algorithm.h"
#include "data_structure/parallel/atomics.h"
#include "data_structure/parallel/random.h"
#include "data_structure/parallel/thread_pool.h"

namespace parallel {

template <typename T>
class thread_container {
public:
        using iterator = typename std::vector<T>::iterator;
        using const_iterator = typename std::vector<T>::const_iterator;

        explicit thread_container(size_t size = 100) {
                m_elems.reserve(size);
                m_counter.store(0, std::memory_order_relaxed);
        }

        inline bool try_pop(T& elem) {
                size_t offset = m_counter.fetch_add(1, std::memory_order_relaxed);

                if (offset < m_elems.size()) {
                        elem = m_elems[offset];
                        return true;
                } else {
                        return false;
                }
        }

        inline void reserve(size_t size) {
                m_elems.reserve(size);
        }

        inline void concurrent_push_back(const T& elem) {
                std::lock_guard<spin_lock> guard(m_lock);
                m_elems.push_back(elem);
        }

        inline void push_back(const T& elem) {
                m_elems.push_back(elem);
        }

        inline void clear() {
                m_elems.clear();
                m_counter.store(0, std::memory_order_relaxed);
        }

        inline size_t size() const {
                return m_elems.size();
        }

        inline bool empty() const {
                return m_counter.load(std::memory_order_relaxed) >= m_elems.size();
        }

        inline void swap(thread_container<T>& container) {
                m_elems.swap(container.m_elems);
                std::swap(m_counter, container.m_counter);
        }

        inline iterator begin() {
                return m_elems.begin();
        }

        inline iterator end() {
                return m_elems.end();
        }

        inline const_iterator begin() const {
                return m_elems.begin();
        }

        inline const_iterator end() const {
                return m_elems.end();
        }

        template <typename... arg_type>
        void emplace_back(arg_type&&... args) {
                m_elems.emplace_back(std::forward<arg_type>(args)...);
        }

        template <typename... arg_type>
        void concurrent_emplace_back(arg_type&&... args) {
                std::lock_guard<spin_lock> guard(m_lock);
                m_elems.emplace_back(std::forward<arg_type>(args)...);
        }
private:
        std::vector<T> m_elems;
        double dummy[parallel::g_cache_line_size - sizeof(std::vector<T>)];
        AtomicWrapper<size_t> m_counter;
        spin_lock m_lock;
};

template <typename T>
class task_queue {
public:
        using iterator = typename Cvector<thread_container<T>>::iterator;

        explicit task_queue(size_t size)
                :       m_thread_containers(size)
        {
                m_counter.store(0, std::memory_order_relaxed);
        }

        inline bool try_pop(T& elem) {
                size_t offset = m_counter.load(std::memory_order_relaxed);

                bool res = false;
                while (true) {
                        offset = m_counter.load(std::memory_order_relaxed);

                        if (offset >= m_thread_containers.size() ||
                                (res = m_thread_containers[offset].get().try_pop(elem))) {
                                break;
                        }

                        m_counter.compare_exchange_strong(offset, offset + 1, std::memory_order_release);
                }
                return res;
        }

        inline bool empty() const {
                return m_counter.load(std::memory_order_relaxed) >= m_thread_containers.size()
                       || is_size_zero();
        }

        inline thread_container<T>& operator[](size_t thread_id) {
                ALWAYS_ASSERT(thread_id < m_thread_containers.size());
                return m_thread_containers[thread_id].get();
        }

        inline const thread_container<T>& operator[](size_t thread_id) const {
                ALWAYS_ASSERT(thread_id < m_thread_containers.size());
                return m_thread_containers[thread_id].get();
        }

        inline void push(const T& elem) {
                m_thread_containers[0].get().push_back(elem);
        }

        inline void clear() {
                for (auto& container : m_thread_containers) {
                        container.get().clear();
                }
                m_counter.store(0, std::memory_order_relaxed);
        }

        inline size_t size() const {
                size_t size = 0;
                for (const auto& container : m_thread_containers) {
                        size += container.get().size();
                }
                return size;
        }

        inline iterator begin() {
                return m_thread_containers.begin();
        }

        inline iterator end() {
                return m_thread_containers.end();
        }

protected:

        inline bool is_size_zero() const {
                for (const auto& container : m_thread_containers) {
                        if (container.get().size() > 0) {
                                return false;
                        }
                }
                return true;
        }

        Cvector<thread_container<T>> m_thread_containers;
        double dummy[parallel::g_cache_line_size - sizeof(Cvector<thread_container<T>>)];
        std::atomic<size_t> m_counter;
};

template <typename T>
class circular_task_queue {
public:
        using iterator = typename Cvector<thread_container<T>>::iterator;

        explicit circular_task_queue(size_t size)
                :       m_thread_containers(size)
        {}

        inline bool try_pop(T& elem, uint32_t id) {
                if (m_thread_containers[id].get().try_pop(elem)) {
                        return true;
                }

                for (uint32_t next_id = id + 1; next_id < m_thread_containers.size(); ++next_id) {
                        if (m_thread_containers[next_id].get().try_pop(elem)) {
                                return true;
                        }
                }

                for (uint32_t next_id = 0; next_id < id; ++next_id) {
                        if (m_thread_containers[next_id].get().try_pop(elem)) {
                                return true;
                        }
                }

                return false;
        }

        inline bool empty() const {
                return is_size_zero();
        }

        inline thread_container<T>& operator[](size_t thread_id) {
                ALWAYS_ASSERT(thread_id < m_thread_containers.size());
                return m_thread_containers[thread_id].get();
        }

        inline const thread_container<T>& operator[](size_t thread_id) const {
                ALWAYS_ASSERT(thread_id < m_thread_containers.size());
                return m_thread_containers[thread_id].get();
        }

        inline void push(const T& elem) {
                m_thread_containers[0].get().push_back(elem);
        }

        inline void clear() {
                for (auto& container : m_thread_containers) {
                        container.get().clear();
                }
        }

        inline size_t size() const {
                size_t size = 0;
                for (const auto& container : m_thread_containers) {
                        size += container.get().size();
                }
                return size;
        }

        inline iterator begin() {
                return m_thread_containers.begin();
        }

        inline iterator end() {
                return m_thread_containers.end();
        }

private:

        inline bool is_size_zero() const {
                for (const auto& container : m_thread_containers) {
                        if (container.get().size() > 0) {
                                return false;
                        }
                }
                return true;
        }

        Cvector<thread_container<T>> m_thread_containers;
        double dummy[parallel::g_cache_line_size - sizeof(Cvector<thread_container<T>>)];
};

static void test_task_queue(size_t num_tests) {
        random rnd(256);
        uint32_t num_threads = g_thread_pool.NumThreads() + 1;

        std::vector<random> thread_rnds;
        thread_rnds.reserve(num_threads);
        for (size_t id = 0; id < num_threads; ++id) {
                thread_rnds.emplace_back(id);
        }
        task_queue<int> queue(num_threads);

        for (size_t num_test = 0; num_test < num_tests; ++num_test) {
                ALWAYS_ASSERT(queue.empty());

                std::cout << "Test " << num_test << " task_queue with " << num_threads << " threads" << std::endl;

                std::vector<std::vector<int>> input(num_threads);
                size_t size_per_thread_min = 1000;
                size_t size_per_thread_max = 5000;

                // generate data
                for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
                        size_t size = rnd.random_number(size_per_thread_min, size_per_thread_max);
                        input[thread_id].reserve(size);
                        for (size_t i = 0; i < size; ++i) {
                                input[thread_id].push_back(rnd.random_number(0, 100000000));
                        }
                }

                // init queue
                std::atomic<uint32_t> counter(0);
                auto insert_task = [&](uint32_t thread_id) {
                        for (auto elem : input[thread_id]) {
                                queue[thread_id].push_back(elem);
                        }
                        thread_rnds[thread_id].shuffle(queue[thread_id].begin(), queue[thread_id].end());
                };

                std::vector<std::future<void>> futures;
                futures.reserve(num_threads);
                for (size_t id = 1; id < num_threads; ++id) {
                        futures.push_back(g_thread_pool.Submit(id - 1, insert_task, id));
                        //futures.push_back(g_thread_pool.Submit(insert_task, id));
                }
                insert_task(0);

                std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                        future.get();
                });
                rnd.shuffle(queue.begin(), queue.end());
                ALWAYS_ASSERT(!queue.empty());

                // read data from task queue

                std::vector<std::vector<int>> output(num_threads);

                auto read_task = [&](uint32_t thread_id) {
                        int elem;
                        while (queue.try_pop(elem)) {
                              output[thread_id].push_back(elem);
                        }
                };
                futures.clear();
                for (size_t id = 1; id < num_threads; ++id) {
                        futures.push_back(g_thread_pool.Submit(id - 1, read_task, id));
                        //futures.push_back(g_thread_pool.Submit(read_task, id));
                }
                read_task(0);

                std::for_each(futures.begin(), futures.end(), [&](auto& future) {
                        future.get();
                });

                std::vector<int> input_overall, output_overall;
                for (size_t i = 0; i < num_threads; ++i) {
                        std::copy(input[i].begin(), input[i].end(), std::back_inserter(input_overall));
                        std::copy(output[i].begin(), output[i].end(), std::back_inserter(output_overall));
                }
                std::sort(input_overall.begin(), input_overall.end());
                std::sort(output_overall.begin(), output_overall.end());

                ALWAYS_ASSERT(input_overall == output_overall);

                queue.clear();
        }

}

};

namespace std {
template <typename T>
void swap(parallel::thread_container<T>& lhs, parallel::thread_container<T>& rhs) {
        lhs.swap(rhs);
}

};