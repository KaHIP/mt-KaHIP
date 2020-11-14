mt-KaHIP v1.00  [![Codacy Badge](https://api.codacy.com/project/badge/Grade/b4dc07be37694dc4841beba6ba038c19)](https://app.codacy.com/manual/schulzchristian/KaHIP?utm_source=github.com&utm_medium=referral&utm_content=schulzchristian/KaHIP&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/KaHIP/KaHIP.svg?branch=master)](https://travis-ci.org//KaHIP) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FOSSA Status]([![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FKaHIP%2Fmt-KaHIP.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FKaHIP%2Fmt-KaHIP?ref=badge_shield))
=====

The graph partitioning framework mt-KaHIP -- Multi-Threaded Karlsruhe High Quality Partitioning.

The graph partitioning problem asks for a division of a graph's node set into k equally sized blocks such that the number of edges that run between the blocks is minimized. KaHIP is a family of graph partitioning programs. Here, we separately release the shared-memory parallel multilevel version of KaHIP which is due to Yaroslav Akhremtsev. Here is an overview of KaHIP:

<p align="center">
<img src="./img/MGPall_en_new.png"
  alt="framework overview"
  width="601" height="558">
</p>

## Main project site:
https://kahip.github.io

Installation Notes
=====
## Downloading mt-KaHIP: 
You can download mt-KaHIP with the following command line:

```console
git clone https://github.com/KaHIP/mt-KaHIP
```

## Compiling mt-KaHIP: 
Before you can start you need to install the following software packages:
- TBB (https://www.threadingbuildingblocks.org/) 
- OpenMP
- CMake

Once you installed the packages, just type 
```console
./compile_withcmake.sh. 
```
In this case, all binaries, libraries and headers are in the folder ./deploy/ 

Alternatively use the standard cmake build process:
```console 
mkdir build
cd build 
cmake ../ -DCMAKE_BUILD_TYPE=Release     
make 
cd ..
```
In this case, the binaries, libraries and headers are in the folder ./build

Running Programs
=====

For a description of the graph format (and an extensive description of all other programs) please have a look into the manual. We give a short examples here.

## Overview of Programs and Usecase

### Default Partitioning Problem  
These programs and configurations take a graph and partition it more or less sequentially. We list here kaffpa and kaffpaE (the evolutionary framework) and their configurations. In general, the configurations are such that you can invest a lot of time into solution quality using the memetic algorithm. The memetic algorithm can also be run in parallel using MPI. In general, the more time and ressources you invest, the better will be the quality of your partition. We have a lot of trade-offs, contact us if you are unsure what works best for your application. For a description of the algorithm have a look at the references that we list in the manual.

| Use Case | Input |  Programs |
| ------------ | -------- | -------- |
| Graph Format || graphchecker  |
| Evaluate Partitions || evaluator |
| Fast Partitioning |  | mtkahip preconfiguration set to socialparallel  |
| Good Partitioning |  | mtkahip preconfiguration set to multitrysocialparallel  |

#### Example Runs
```console
./deploy/graphchecker ./examples/rgg_n_2_15_s0.graph 
```

```console
./deploy/mtkahip ./examples/rgg_n_2_15_s0.graph --k 4  --preconfiguration=socialparallel --num_threads=4
```

Linking the mt-KaHIP Library 
=====
mt-KaHIP also offers libaries and interfaces to link the algorithms directly to your code. We explain the details of the interface in the manual. Below we list an example program that links the kahip library. This example can also  be found in misc/example_library_call/. The program will also be compiled when running the cmake process from above (see CMakeLists.txt file for necessary libraries etc...). 

```cpp
#include <iostream>
#include <sstream>

#include "mtkahip_interface.h"


int main(int argn, char **argv) {

        std::cout <<  "partitioning graph from the manual"  << std::endl;

        int n            = 5;
        int* xadj        = new int[6];
        xadj[0] = 0; xadj[1] = 2; xadj[2] = 5; xadj[3] = 7; xadj[4] = 9; xadj[5] = 12;

        int* adjncy      = new int[12];
        adjncy[0]  = 1; adjncy[1]  = 4; adjncy[2]  = 0; adjncy[3]  = 2; adjncy[4]  = 4; adjncy[5]  = 1; 
        adjncy[6]  = 3; adjncy[7]  = 2; adjncy[8]  = 4; adjncy[9]  = 0; adjncy[10] = 1; adjncy[11] = 3; 
        
        double imbalance = 0.03;
        int* part        = new int[n];
        int edge_cut     = 0;
        int nparts       = 2;
        int numthreads   = 4;
        int* vwgt        = NULL;
        int* adjcwgt     = NULL;

        mtkahip(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, false, 0, FASTSOCIALMULTITRY_PARALLEL, numthreads, & edge_cut, part);

        std::cout <<  "edge cut " <<  edge_cut  << std::endl;
                
}
```

Licence
=====
The program is licenced under MIT licence.
If you publish results using our algorithms, please acknowledge our work by quoting the following paper:

```
@article{DBLP:journals/tpds/AkhremtsevSS20,
  author    = {Yaroslav Akhremtsev and
               Peter Sanders and
               Christian Schulz},
  title     = {High-Quality Shared-Memory Graph Partitioning},
  journal   = {{IEEE} Trans. Parallel Distributed Syst.},
  volume    = {31},
  number    = {11},
  pages     = {2710--2722},
  year      = {2020},
  url       = {https://doi.org/10.1109/TPDS.2020.3001645},
  doi       = {10.1109/TPDS.2020.3001645},
  timestamp = {Fri, 02 Oct 2020 14:40:05 +0200},
  biburl    = {https://dblp.org/rec/journals/tpds/AkhremtsevSS20.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

mt-KaHIP Project Contributors (sorted by last name)
=====
Yaroslav Akhremtsev

Peter Sanders

Christian Schulz
