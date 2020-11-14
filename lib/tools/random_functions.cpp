/******************************************************************************
 * random_functions.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include "random_functions.h"

thread_local MersenneTwister random_functions::m_mt;
thread_local int random_functions::m_seed = 0;

random_functions::random_functions()  {
}

random_functions::~random_functions() {
}
