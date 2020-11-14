#pragma once

#include "data_structure/parallel/atomics.h"

namespace parallel {

class spin_lock {

public:
        spin_lock()
                : m_locked(false) {}

        spin_lock(const spin_lock&) = default;
        spin_lock(spin_lock&&) = default;

        spin_lock& operator=(const spin_lock&) = default;
        spin_lock& operator=(spin_lock&&) = default;

        inline void lock() {
                bool expected = false;
                while (!m_locked.compare_exchange_weak(expected, true, std::memory_order_acquire)) {
                        expected = false;
                }
        }

        inline bool try_lock() {
                bool expected = false;
                return m_locked.compare_exchange_strong(expected, true, std::memory_order_acquire);
        }

        inline void unlock() {
                m_locked.store(false, std::memory_order_release);
        }

private:
        parallel::AtomicWrapper<bool> m_locked;
};

}