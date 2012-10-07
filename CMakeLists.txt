# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# TODO: Add install code.

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

project(octopus CXX)

# Project version.
set(OCTOPUS_MAJOR_VERSION 0)
set(OCTOPUS_MINOR_VERSION 1)
set(OCTOPUS_PATCH_LEVEL   0)
set(OCTOPUS_SOVERSION     ${OCTOPUS_MAJOR_VERSION})
set(OCTOPUS_VERSION
    "${OCTOPUS_MAJOR_VERSION}.${OCTOPUS_MINOR_VERSION}.${OCTOPUS_PATCH_LEVEL}")

set(CMAKE_MODULE_PATH ${octopus_SOURCE_DIR}/cmake
    ${HPX_ROOT}/share/cmake-${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}/Modules)

if(MSVC)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

find_package(HPX REQUIRED)

if(NOT MSVC)
  set(HPX_RPATH
      "${HPX_RPATH}:${CMAKE_INSTALL_PREFIX}/lib/hpx:${CMAKE_BINARY_DIR}/lib/hpx")
endif()

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})

include_directories(${octopus_SOURCE_DIR})

################################################################################
find_package(HPX_Silo REQUIRED)

include_directories(${SILO_INCLUDE_DIR})
link_directories(${SILO_LIBRARY_DIR})

################################################################################
hpx_option(OCTOPUS_VERIFY BOOL
           "Enable verification of Octopus (default: ON)" ON ADVANCED)

if(OCTOPUS_VERIFY)
  add_definitions(-DOCTOPUS_ENABLE_VERIFICATION=1)
else()
  add_definitions(-DOCTOPUS_ENABLE_VERIFICATION=0)
endif()

################################################################################
hpx_option(OCTOPUS_TEST_IN_PLACE BOOL
           "Enable in-code testing of Octopus (default: ON)" ON ADVANCED)

if(OCTOPUS_TEST_IN_PLACE)
  add_definitions(-DOCTOPUS_ENABLE_TEST_IN_PLACE=1)
else()
  add_definitions(-DOCTOPUS_ENABLE_TEST_IN_PLACE=0)
endif()

################################################################################
add_subdirectory(src)

################################################################################
add_hpx_pseudo_target(tests)
add_subdirectory(tests)

################################################################################
add_hpx_pseudo_target(applications)
add_subdirectory(applications)


