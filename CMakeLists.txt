cmake_minimum_required(VERSION 3.16)
project(KaHIP LANGUAGES C CXX)

include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache: ${CCACHE_PROGRAM}")
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# If no build mode is specified, build in Release mode.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(64BITMODE "64 bit mode" OFF)
option(OPTIMIZED_OUTPUT "optimized output" OFF)

function(add_supported_cxx_option option_name)
  string(MAKE_C_IDENTIFIER "${option_name}" option_id)
  set(test_var "COMPILER_SUPPORTS_${option_id}")
  check_cxx_compiler_flag("${option_name}" "${test_var}")
  if(${test_var})
    target_compile_options(mtkahip_obj PRIVATE "${option_name}")
  endif()
endfunction()

# Check dependencies
find_package(OpenMP REQUIRED)
find_package(TBB QUIET COMPONENTS tbb tbbmalloc tbbmalloc_proxy)
find_package(Threads REQUIRED)

find_path(NUMA_INCLUDE_DIR NAMES numa.h)
find_library(NUMA_LIBRARY NAMES numa REQUIRED)

if(TBB_FOUND)
  set(MTKAHIP_TBB_LIBS TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)
else()
  set(MTKAHIP_TBB_LIBS tbb tbbmalloc tbbmalloc_proxy)
endif()

set(MTKAHIP_COMMON_LIBS
  OpenMP::OpenMP_CXX
  ${NUMA_LIBRARY}
  ${MTKAHIP_TBB_LIBS}
  Threads::Threads
  atomic
  dl
)

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
  extern/argtable3-3.0.3/argtable3.c
)

add_library(mtkahip_obj OBJECT ${LIBKAFFPA_SOURCE_FILES})
target_include_directories(mtkahip_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/app
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/io
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/partition
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/tools
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/growt
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/ips4o
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/stats/include
    $<$<BOOL:${NUMA_INCLUDE_DIR}>:${NUMA_INCLUDE_DIR}>
)
target_compile_definitions(mtkahip_obj
  PUBLIC
    CPP11THREADS
    _REENTRANT
    MODE_KAFFPA
    POINTER64=1
    GIT_DESC="942944c"
)

if(64BITMODE)
  target_compile_definitions(mtkahip_obj PUBLIC MODE64BITEDGES POINTER64=1)
endif()

if(OPTIMIZED_OUTPUT)
  target_compile_definitions(mtkahip_obj PUBLIC KAFFPAOUTPUT)
endif()

add_supported_cxx_option(-funroll-loops)
add_supported_cxx_option(-mcx16)
add_supported_cxx_option(-Wno-unused-but-set-variable)
add_supported_cxx_option(-fno-stack-limit)
add_supported_cxx_option(-Wall)
add_supported_cxx_option(-march=native)
add_supported_cxx_option(-fpermissive)
add_supported_cxx_option(-Wno-unused-result)
add_supported_cxx_option(-Wno-sign-compare)
add_supported_cxx_option(-Wno-unused-function)

# generate targets for each binary
add_executable(mtkahip app/mtkahip.cpp)
target_compile_definitions(mtkahip PRIVATE MODE_KAFFPA)
target_link_libraries(mtkahip PRIVATE mtkahip_obj ${MTKAHIP_COMMON_LIBS})
install(TARGETS mtkahip DESTINATION bin)

# Shared interface library
add_library(interfacemtkahip SHARED interface/mtkahip_interface.cpp)
target_include_directories(interfacemtkahip PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interfacemtkahip PRIVATE MODE_KAFFPA)
target_link_libraries(interfacemtkahip PUBLIC mtkahip_obj ${MTKAHIP_COMMON_LIBS})
if(LIB_METIS)
  target_link_libraries(interfacemtkahip PUBLIC ${LIB_METIS})
endif()
install(TARGETS interfacemtkahip DESTINATION lib)

# Static interface library
add_library(interfacemtkahip_static STATIC interface/mtkahip_interface.cpp)
target_include_directories(interfacemtkahip_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interfacemtkahip_static PRIVATE MODE_KAFFPA)
target_link_libraries(interfacemtkahip_static PUBLIC mtkahip_obj ${MTKAHIP_COMMON_LIBS})
set_target_properties(interfacemtkahip_static PROPERTIES PUBLIC_HEADER interface/mtkahip_interface.h)

add_executable(graphchecker app/graphchecker.cpp)
target_compile_definitions(graphchecker PRIVATE MODE_GRAPHCHECKER)
target_link_libraries(graphchecker PRIVATE mtkahip_obj ${MTKAHIP_COMMON_LIBS})
install(TARGETS graphchecker DESTINATION bin)

add_executable(evaluator app/evaluator.cpp)
target_compile_definitions(evaluator PRIVATE MODE_EVALUATOR)
target_link_libraries(evaluator PRIVATE mtkahip_obj ${MTKAHIP_COMMON_LIBS})
install(TARGETS evaluator DESTINATION bin)

add_executable(interface_test misc/example_library_call/interface_test.cpp interface/mtkahip_interface.cpp)
target_include_directories(interface_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interface_test PRIVATE MODE_KAFFPA)
target_link_libraries(interface_test PRIVATE mtkahip_obj ${MTKAHIP_COMMON_LIBS})
if(LIB_METIS)
  target_link_libraries(interface_test PRIVATE ${LIB_METIS})
endif()
install(TARGETS interface_test DESTINATION bin)

install(TARGETS interfacemtkahip_static
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
