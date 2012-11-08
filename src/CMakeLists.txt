# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

add_definitions(-DOCTOPUS_INTERNAL)
add_definitions(-DOCTOPUS_EXPORTS)

add_hpx_component(octopus
  SOURCES
          octopus_component.cpp
          driver.cpp
          default_main.cpp
          child_index.cpp
          engine/engine_interface.cpp
          engine/engine_server.cpp
          engine/runtime_config.cpp
          octree/octree_client.cpp
          octree/octree_server.cpp
          science/reconstruction.cpp
          science/science_table.cpp
          io/silo.cpp
  DEPENDENCIES ${SILO_LIBRARY}
  FOLDER "Core"
  ESSENTIAL)

if(VISIT_DRIVER_FOUND AND VISIT_LIBSIM_FOUND)
  add_hpx_component(visit_simulation
    SOURCES
            visit/visit_component.cpp
            visit/visit_simulation_server.cpp
    DEPENDENCIES ${VISIT_LIBSIM_LIBRARY}
    COMPONENT_DEPENDENCIES octopus
    FOLDER "Core"
    ESSENTIAL)
endif()

