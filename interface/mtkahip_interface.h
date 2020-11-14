/******************************************************************************
 * mtkahip_interface.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
