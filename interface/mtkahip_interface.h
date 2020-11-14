/******************************************************************************
 * mtkahip_interface.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/


#ifndef MTKAHIP_INTERFACE_RYEEZ6WJ
#define MTKAHIP_INTERFACE_RYEEZ6WJ

#include <cstdint>

#ifdef __cplusplus

extern "C"
{
#endif

const int FASTSOCIALMULTITRY_PARALLEL = 6;
const int FASTSOCIAL_PARALLEL = 7;

// same data structures as in metis 
// edgecut and part are output parameters
// part has to be an array of n ints
void mtkahip(int* n, int* vwgt, int* xadj, 
                   int* adjcwgt, int* adjncy, int* nparts, 
                   double* imbalance,  bool suppress_output, int seed, int mode, uint32_t num_threads,
                   int* edgecut, int* part);

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: KAFFPA_INTERFACE_RYEEZ6WJ */
