#pragma once

#include <cstdint>

namespace parallel {

inline static uint8_t most_significant_bit_index(size_t num) {
        uint8_t res = 0;

        while (num >>= 1) {
                ++res;
        }

        return res;
}

inline static uint8_t log2(size_t num) {
        return most_significant_bit_index(num);
}

inline static uint8_t bits_number(size_t num) {
        return most_significant_bit_index(num) + 1;
}

constexpr static uint32_t round_up_to_next_power_2(uint32_t v) {
        v--;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;
        return v;
}

constexpr static uint32_t round_up_to_previous_power_2(uint32_t x) {
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >> 16);
        return x - (x >> 1);
}

constexpr static bool is_power_2(uint64_t x) {
        return ((x != 0) && !(x & (x - 1)));
}
}
