#pragma once

#include <atomic>
#include <iostream>

namespace parallel {

template<typename T>
struct AtomicWrapper;

template<typename T>
std::ostream& operator<<(std::ostream& out, const AtomicWrapper<T>& atomic);

template<typename T>
struct AtomicWrapper {
        std::atomic <T> m_atomic;

        inline AtomicWrapper(T value = T())
                : m_atomic(value) {}

        inline explicit AtomicWrapper(const std::atomic <T>& atomic)
                : m_atomic(atomic.load(std::memory_order_acquire)) {}

        inline explicit AtomicWrapper(const AtomicWrapper& other)
                : m_atomic(other.m_atomic.load(std::memory_order_acquire)) {}

        inline AtomicWrapper& operator=(const AtomicWrapper& other) {
                m_atomic.store(other.m_atomic.load(std::memory_order_acquire), std::memory_order_release);
                return *this;
        }

        inline AtomicWrapper& operator=(const T& other) {
                m_atomic.store(other);
                return *this;
        }

        inline AtomicWrapper& operator+= (T other) {
                m_atomic.fetch_add(other);
                return *this;
        }

        inline bool operator!() const {
                return !m_atomic;
        }

        inline T load(std::memory_order memory_order = std::memory_order_seq_cst) const {
                return m_atomic.load(memory_order);
        }

        inline void store(const T& other, std::memory_order memory_order = std::memory_order_seq_cst) {
                m_atomic.store(other, memory_order);
        }

        inline AtomicWrapper& operator++() {
                m_atomic.fetch_add(1, std::memory_order_acq_rel);
                return *this;
        }

        inline bool operator<(const T& other) const {
                return m_atomic < other;
        }

        inline bool operator<(const AtomicWrapper& other) const {
                return m_atomic < other.m_atomic;
        }

        inline bool operator>(const T& other) const {
                return m_atomic > other;
        }

        inline bool operator>(const AtomicWrapper& other) const {
                return m_atomic > other.m_atomic;
        }

        inline bool operator==(const T& value) const {
                return m_atomic == value;
        }

        inline operator T() const {
                return m_atomic.load(std::memory_order_seq_cst);
        }

        inline T fetch_add(T value, std::memory_order order = std::memory_order_seq_cst) {
                return m_atomic.fetch_add(value, order);
        }

        inline T fetch_sub(T value, std::memory_order order = std::memory_order_seq_cst) {
                return m_atomic.fetch_sub(value, order);
        }

        inline T exchange(T desired, std::memory_order order = std::memory_order_seq_cst) {
                return m_atomic.exchange(desired, order);
        }

        inline bool compare_exchange_weak(T& expected, T desired,
                                          std::memory_order memory_order = std::memory_order_seq_cst) {
                return m_atomic.compare_exchange_weak(expected, desired, memory_order);
        }

        inline bool
        compare_exchange_weak(T& expected, T desired, std::memory_order success, std::memory_order failure) {
                return m_atomic.compare_exchange_weak(expected, desired, success, failure);
        }

        inline bool compare_exchange_strong(T& expected, T desired,
                                            std::memory_order memory_order = std::memory_order_seq_cst) {
                return m_atomic.compare_exchange_strong(expected, desired, memory_order);
        }

        inline bool
        compare_exchange_strong(T& expected, T desired, std::memory_order success, std::memory_order failure) {
                return m_atomic.compare_exchange_strong(expected, desired, success, failure);
        }

        friend std::ostream& operator<<<>(std::ostream& out, const AtomicWrapper<T>& atomic);
};

template<typename T>
std::ostream& operator<<(std::ostream& out, const AtomicWrapper<T>& atomic) {
        return out << atomic.m_atomic.load(std::memory_order_acquire);
}

}