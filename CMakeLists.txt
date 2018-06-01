#########################
# Important note: This CMake build file is intended for using MiniSat
# as a static library. It does not build the DIMACS-parsing binaries,
# and is only intended for consumption by other CMake projects.
#########################
cmake_minimum_required(VERSION 3.5)
project(MiniSat VERSION 2.2 LANGUAGES CXX)

include(GNUInstallDirs)


add_library(libminisat STATIC
    minisat/core/Solver.cc
    minisat/utils/Options.cc
    minisat/utils/System.cc
    minisat/simp/SimpSolver.cc
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

# ToDo: With some modernization of MiniSat, this is not strictly neccessary
target_compile_definitions(libminisat PUBLIC __STDC_LIMIT_MACROS __STDC_FORMAT_MACROS)


# Also build two MiniSat executables
add_executable(minisat
    minisat/core/Main.cc
)
target_link_libraries(minisat libminisat)


add_executable(minisat-simp
    minisat/simp/Main.cc
)
target_link_libraries(minisat-simp libminisat)

# Workaround for libstdc++ + Clang + -std=gnu++11 bug.
set_target_properties(libminisat minisat minisat-simp
    PROPERTIES
      CXX_EXTENSIONS OFF
)


if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
    target_compile_options( libminisat PRIVATE -Wall -Wextra )
    target_compile_options( minisat PRIVATE -Wall -Wextra )
    target_compile_options( minisat-simp PRIVATE -Wall -Wextra )
endif()
if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    target_compile_options( libminisat PRIVATE /W4 )
    target_compile_options( minisat PRIVATE /W4 )
    target_compile_options( minisat-simp PRIVATE /W4 )
endif()


###############
# Installation
##
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/MiniSat)

install(
    TARGETS
      libminisat
    EXPORT MiniSatTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    TARGETS
      minisat
      minisat-simp
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


install(DIRECTORY minisat/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/minisat FILES_MATCHING PATTERN "*.h*")

install(EXPORT MiniSatTargets
    FILE MiniSatTargets.cmake
    NAMESPACE MiniSat::
    DESTINATION ${INSTALL_CONFIGDIR}
)


#####################
# ConfigVersion file
##
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/CMake/MiniSatConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

## Install all the helper files
install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfig.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/MiniSatConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
