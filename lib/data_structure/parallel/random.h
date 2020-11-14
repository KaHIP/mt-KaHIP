#pragma once

#include <random>
#include <iostream>
namespace parallel {
class random {
public:
        explicit random(uint32_t seed)
                :       m_rnd_bit(0, 1)
                ,       m_mt(seed)
        {}

        inline bool bit() {
                return m_rnd_bit(m_mt);
        }

        // returns random number in [a, b]
        template <typename T = uint32_t>
        inline T random_number(T a = std::numeric_limits<T>::min(), T b = std::numeric_limits<T>::max()) {

                using dist_type = typename std::conditional<std::is_integral<T>::value,
                                                            std::uniform_int_distribution<T>,
                                                            std::uniform_real_distribution<T>>::type;
                dist_type rnd(a, b);
                return rnd(m_mt);
        }

        void set_seed(uint32_t seed) {
                m_mt.seed(seed);
        }

        template <typename T>
        inline void shuffle(std::vector<T>& vec) {
                shuffle(vec.begin(), vec.end());
        }

        template <typename iterator_type>
        inline void shuffle(iterator_type begin, iterator_type end) {
                size_t size = end - begin;
                if (size < 2) {
                        return;
                }

                std::uniform_int_distribution<uint32_t> rnd;
                for (size_t i = 0; i < size; ++i) {
                        size_t rnd_ind = i + rnd(m_mt) % (size - i);
                        std::swap(*(begin + i), *(begin + rnd_ind));
                }
        }

        template <typename iterator_type>
        inline void shuffle_blocks(iterator_type begin, iterator_type end, size_t block_size = 4) {
                size_t size = end - begin;
                if(size < 10 || size < block_size) {
                        shuffle(begin, end);
                        return;
                }

                const size_t distance = 100;
                std::uniform_int_distribution<uint32_t> rnd(0, distance);
                for( unsigned int i = 0; i + block_size < size; i++) {
                        size_t pos_a = i;
                        size_t pos_b = (pos_a + rnd(m_mt)) % (size - block_size + 1);

                        for (size_t j = 0; j < block_size; ++j) {
                                std::swap(*(begin + pos_a + j), *(begin + pos_b + j));
                        }
                }
        }
private:
        std::uniform_int_distribution<uint32_t> m_rnd_bit;
        std::mt19937 m_mt;
};
}