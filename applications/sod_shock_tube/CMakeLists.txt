# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

add_definitions(-DOCTOPUS_CURRENT_SOURCE_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}")

set(applications
    sod_shock_tube_driver
   )

foreach(application ${applications})
  set(sources ${application}.cpp)

  source_group("Source Files" FILES ${sources})

  add_hpx_executable(${application} SOURCES ${sources} ${${application}_FLAGS}
                     DEPENDENCIES octopus_component simV2
                     FOLDER "Applications")

  add_hpx_pseudo_target(applications.${application})
  add_hpx_pseudo_dependencies(applications applications.${application})
  add_hpx_pseudo_dependencies(applications.${application} ${application}_exe)
endforeach()

