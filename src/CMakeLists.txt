# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

add_definitions("-DOCTOPUS_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")

add_hpx_component(octopus
  SOURCES component_module.cpp
          driver.cpp
          octree/octree_client.cpp
          octree/octree_server.cpp  
          engine/engine_server.cpp
  COMPONENT_DEPENDENCIES distributing_factory
  FOLDER "Core"
  ESSENTIAL)

add_subdirectory(modules)

