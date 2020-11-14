KaHIP v1.00  
=====  
  
Parallel graph partitioning framework KaHIP -- Multi-Threaded Karlsruhe High Quality Partitioning.  
  
The graph partitioning problem asks for a division of a graph's node set into k equally sized blocks such that the number of edges that run between the blocks is minimized. KaHIP is a family of graph partitioning programs. It includes KaFFPa (Karlsruhe Fast Flow Partitioner), which is a multilevel graph partitioning algorithm, in its variants Strong, Eco and Fast, KaFFPaE (KaFFPaEvolutionary) which is a parallel evolutionary algorithm that uses KaFFPa to provide combine and mutation operations, as well as KaBaPE which extends the evolutionary algorithm. Moreover, specialized techniques are included to partition road networks (Buffoon), to output a vertex separator from a given partition as well as techniques geared towards the efficient partitioning of social networks.  
  
Main project site:  
http://algo2.iti.kit.edu/documents/kahip/index.html  
  
Installation Notes  
=====  
  
Before you can start you need to install the following software packages:  
  
- Scons (http://www.scons.org/)  
- Argtable (http://argtable.sourceforge.net/)
- TBB (https://www.threadingbuildingblocks.org/) 
  
Getting the code:  

    git clone git@github.com:yarchi/KaHIP.git  
    git checkout add_parallel_local_search  
    git submodule init  
    git submodule update  

  
Compiling binary:  

    scons program=kaffpa variant=optimized -j 16  

  
Running the partitioner:  

    ./optimized/kaffpa <graph_name> --k <number of blocks> --preconfiguration=fastsocialmultitry_parallel --num_threads=16 --seed=0 --output_filename graph.part  

  
  
Compiling the library:  

    scons program=library variant=optimized -j 16  

  
Using the library:  
To use the library, you need the header file *interface/kaHIP_interface.h* that contains the declaration of  the function kaffpa.

Compile you binary with the following additional flags: `-ltbb -ltbb_malloc -ltbbmalloc_proxy -lpthread`