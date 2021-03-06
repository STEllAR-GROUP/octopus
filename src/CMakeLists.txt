# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

add_definitions(-DOCTOPUS_INTERNAL)
add_definitions(-DOCTOPUS_EXPORTS)

if(SILO_FOUND)
  add_hpx_component(octopus
    SOURCES
            octopus_component.cpp
            driver.cpp
            child_index.cpp
            engine/engine_interface.cpp
            engine/engine_server.cpp
            engine/runtime_config.cpp
            octree/octree_client.cpp
            octree/octree_server.cpp
            science/minmod_reconstruction.cpp
            science/ppm_reconstruction.cpp
            science/science_table.cpp
            io/silo.cpp
            io/fstream.cpp
    DEPENDENCIES ${SILO_LIBRARY} dl
    FOLDER "Core"
    ESSENTIAL)
else()
  add_hpx_component(octopus
    SOURCES
            octopus_component.cpp
            driver.cpp
            child_index.cpp
            engine/engine_interface.cpp
            engine/engine_server.cpp
            engine/runtime_config.cpp
            octree/octree_client.cpp
            octree/octree_server.cpp
            science/minmod_reconstruction.cpp
            science/ppm_reconstruction.cpp
            science/science_table.cpp
            io/fstream.cpp
    FOLDER "Core"
    ESSENTIAL)
endif()

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

