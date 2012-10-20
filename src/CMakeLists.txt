# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

add_definitions(-DOCTOPUS_INTERNAL)
add_definitions(-DOCTOPUS_EXPORTS)

add_hpx_component(octopus
  SOURCES
          engine/engine_server.cpp
          engine/runtime_config.cpp
          octree/octree_client.cpp
          octree/octree_server.cpp
          component_module.cpp
          science/reconstruction.cpp
          science/science_table.cpp
          driver.cpp
          default_main.cpp
          child_index.cpp
  FOLDER "Core"
  ESSENTIAL)

