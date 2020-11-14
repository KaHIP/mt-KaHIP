#pragma once

#include "data_structure/parallel/random.h"

namespace parallel {
struct thread_config {
public:
        uint32_t id;
        parallel::random rnd;

        thread_config(uint32_t _id,
                      uint32_t _seed)
                : id(_id),
                  rnd(_seed)
        {}

};
}