#########################
# Important note: Minisat does not properly annotate its API with
# dllimport/dllexport, so the library build is forced to be static.
#########################
cmake_minimum_required(VERSION 3.5)
project(MiniSat VERSION 2.2 LANGUAGES CXX)

option(MINISAT_BUILD_TESTING "Build and run MiniSat's tests" ON)
option(MINISAT_TEST_BENCHMARKS "Register benchmarks with CTest" OFF)

include(GNUInstallDirs)

add_library(libminisat STATIC
    # Impl files
    minisat/core/Solver.cc
    minisat/core/SolverTypes.cc
    # Header files for IDEs
    minisat/core/Solver.h
    minisat/core/SolverTypes.h
    minisat/mtl/Heap.h
)

# Keep the library named as either libminisat.a or minisat.lib
# While having the target's logical name be distinct from minisat (the binary)
set_target_properties(libminisat
    PROPERTIES
      OUTPUT_NAME "minisat"
)


## Add a namespace alias.
# This is useful to abstract over use of the library as installed vs subdirectory build
add_library(MiniSat::libminisat ALIAS libminisat)

target_compile_features(libminisat
    PUBLIC
      cxx_attributes
      cxx_defaulted_functions
      cxx_deleted_functions
      cxx_final
)

target_include_directories(libminisat
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Also build two MiniSat executables
add_executable(minisat
    minisat/core/Main.cc
)
target_link_libraries(minisat libminisat)
