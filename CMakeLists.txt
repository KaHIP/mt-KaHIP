cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

#add_definitions("-funroll-loops  -fno-stack-limit -std=c++1z -mcx16")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-mcx16 COMPILER_SUPPORTS_MCX16)
if(COMPILER_SUPPORTS_MCX16)
  add_definitions(-mcx16)
endif()
#CHECK_CXX_COMPILER_FLAG(-std=c++1z COMPILER_SUPPORTS_CPP1Z)
#if(COMPILER_SUPPORTS_CPP1Z)
  #add_definitions(-std=c++1z)
#endif()

CHECK_CXX_COMPILER_FLAG(-Wno-unused-but-set-variable COMPILER_SUPPORTS_NOUSEDBUTSET)
if(COMPILER_SUPPORTS_NOUSEDBUTSET)
  add_definitions(-Wno-unused-but-set-variable)
endif()

CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-unused-result COMPILER_SUPPORTS_UNUSED)
if(COMPILER_SUPPORTS_UNUSED)
  add_definitions(-Wno-unused-result)
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-sign-compare COMPILER_SUPPORTS_NOSIGNCOMP)
if(COMPILER_SUPPORTS_NOSIGNCOMP)
  add_definitions(-Wno-sign-compare)
endif()


add_definitions(-DCPP11THREADS -D_REENTRANT -U_OPENMP -Wno-unused-function -O3 -DNDEBUG -funroll-loops -std=c++1z -DMODE_KAFFPA -DPOINTER64=1 -DGIT_DESC=\"942944c\" )



# check dependencies
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  add_definitions(${OpenMP_CXX_FLAGS})
endif()
#find_package(TBB REQUIRED)
#if(TBB_CXX_FOUND)
  #add_definitions(${TBB_CXX_FLAGS})
#endif()



# 64 Bit option
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()


# optimized output
option(OPTIMIZED_OUTPUT "optimized output" OFF)
if(OPTIMIZED_OUTPUT)
  add_definitions("-DKAFFPAOUTPUT")
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/growt/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/ips4o)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/stats/include)


set(LIBKAFFPA_SOURCE_FILES
  lib/data_structure/graph_hierarchy.cpp
  lib/algorithms/strongly_connected_components.cpp
  lib/algorithms/topological_sort.cpp
  lib/algorithms/push_relabel.cpp
  lib/io/graph_io.cpp
  lib/tools/quality_metrics.cpp
  lib/tools/random_functions.cpp
  lib/tools/graph_extractor.cpp
  lib/tools/misc.cpp
  lib/tools/partition_snapshooter.cpp
  lib/partition/graph_partitioner.cpp
  lib/partition/w_cycles/wcycle_partitioner.cpp
  lib/partition/coarsening/coarsening.cpp
  lib/partition/coarsening/contraction.cpp
  lib/partition/coarsening/edge_rating/edge_ratings.cpp
  lib/partition/coarsening/matching/matching.cpp
  lib/partition/coarsening/matching/random_matching.cpp
  lib/partition/coarsening/matching/gpa/path.cpp
  lib/partition/coarsening/matching/gpa/gpa_matching.cpp
  lib/partition/coarsening/matching/gpa/path_set.cpp
  lib/partition/coarsening/clustering/node_ordering.cpp
  lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp
  lib/partition/initial_partitioning/initial_partitioning.cpp
  lib/partition/initial_partitioning/initial_partitioner.cpp
  lib/partition/initial_partitioning/initial_partition_bipartition.cpp
  lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp
  lib/partition/initial_partitioning/bipartition.cpp
  lib/partition/initial_partitioning/initial_node_separator.cpp
  lib/partition/uncoarsening/uncoarsening.cpp
  lib/partition/uncoarsening/separator/area_bfs.cpp
  lib/partition/uncoarsening/separator/vertex_separator_algorithm.cpp
  lib/partition/uncoarsening/separator/vertex_separator_flow_solver.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/problem_factory.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph.cpp
  lib/partition/uncoarsening/refinement/mixed_refinement.cpp
  lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp
  lib/partition/uncoarsening/refinement/refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/two_way_flow_refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/boundary_bfs.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/cut_flow_problem_solver.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/most_balanced_minimum_cuts/most_balanced_minimum_cuts.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/quotient_graph_scheduling.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/active_block_quotient_graph_scheduler.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph_fabric.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.cpp
  lib/partition/uncoarsening/refinement/node_separators/greedy_ns_local_search.cpp
  lib/partition/uncoarsening/refinement/node_separators/fm_ns_local_search.cpp
  lib/partition/uncoarsening/refinement/node_separators/localized_fm_ns_local_search.cpp
  lib/algorithms/cycle_search.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.cpp
  lib/partition/uncoarsening/refinement/tabu_search/tabu_search.cpp
  lib/partition/coarsening/matching/local_max.cpp
  lib/data_structure/parallel/thread_pool.cpp
  lib/partition/coarsening/min_hash/hash_common_neighborhood.cpp
  lib/partition/initial_partitioning/parallel/initial_partitioning.cpp
  lib/partition/uncoarsening/refinement/parallel_kway_graph_refinement/multitry_kway_fm.cpp
  lib/partition/uncoarsening/refinement/parallel_kway_graph_refinement/kway_graph_refinement_core.cpp
  lib/partition/uncoarsening/parallel_uncoarsening.cpp
  extern/argtable3-3.0.3/argtable3.c)
add_library(libmtkahip OBJECT ${LIBKAFFPA_SOURCE_FILES})


# generate targets for each binary
add_executable(mtkahip app/mtkahip.cpp $<TARGET_OBJECTS:libmtkahip>)
target_compile_definitions(mtkahip PRIVATE "-DMODE_KAFFPA")
target_link_libraries(mtkahip ${OpenMP_CXX_LIBRARIES} "tbb" "tbbmalloc" "tbbmalloc_proxy" "atomic" "dl" "numa")
install(TARGETS mtkahip DESTINATION bin)

# Shared interface library
add_library(interfacemtkahip SHARED interface/mtkahip_interface.cpp $<TARGET_OBJECTS:libmtkahip>)
target_include_directories(interfacemtkahip PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interfacemtkahip PRIVATE "-DMODE_KAFFPA")
target_link_libraries(interfacemtkahip PUBLIC ${OpenMP_CXX_LIBRARIES} "tbb" "tbbmalloc" "tbbmalloc_proxy" "atomic" "dl" "numa")
if(LIB_METIS)
        target_link_libraries(interface PUBLIC ${LIB_METIS})
endif()
install(TARGETS interfacemtkahip DESTINATION lib)

# Static interface library
add_library(interfacemtkahip_static interface/mtkahip_interface.cpp $<TARGET_OBJECTS:libmtkahip>)
target_include_directories(interfacemtkahip_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interfacemtkahip_static PRIVATE "-DMODE_KAFFPA")
target_link_libraries(interfacemtkahip_static PUBLIC ${OpenMP_CXX_LIBRARIES}  "tbb" "tbbmalloc" "tbbmalloc_proxy" "atomic" "dl" "numa")
set_target_properties(interfacemtkahip_static PROPERTIES PUBLIC_HEADER interface/mtkahip_interface.h)

add_executable(graphchecker app/graphchecker.cpp $<TARGET_OBJECTS:libmtkahip>)
target_compile_definitions(graphchecker PRIVATE "-DMODE_GRAPHCHECKER"  )
target_link_libraries(graphchecker ${OpenMP_CXX_LIBRARIES} "tbb" "tbbmalloc" "tbbmalloc_proxy" "atomic" "dl" "numa")
install(TARGETS graphchecker DESTINATION bin)

add_executable(evaluator app/evaluator.cpp $<TARGET_OBJECTS:libmtkahip>)
target_compile_definitions(evaluator PRIVATE "-DMODE_EVALUATOR")
target_link_libraries(evaluator ${OpenMP_CXX_LIBRARIES}  "tbb" "tbbmalloc" "tbbmalloc_proxy" "atomic" "dl" "numa")
install(TARGETS evaluator DESTINATION bin)

add_executable(interface_test misc/example_library_call/interface_test.cpp interface/mtkahip_interface.cpp $<TARGET_OBJECTS:libmtkahip> )
target_include_directories(interface_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interface_test PRIVATE "-DMODE_KAFFPA")
target_link_libraries(interface_test ${OpenMP_CXX_LIBRARIES} "tbb" "tbbmalloc" "tbbmalloc_proxy" "atomic" "dl" "numa")
if(LIB_METIS)
        target_link_libraries(interface_test ${LIB_METIS})
endif()
install(TARGETS interface_test DESTINATION bin)

install(TARGETS interfacemtkahip_static 
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )



