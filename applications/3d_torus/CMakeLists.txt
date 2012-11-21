# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

add_definitions(-DOCTOPUS_CURRENT_SOURCE_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}")

set(applications
    3d_torus_driver
    3d_torus_update_kappa
    3d_torus_cli
    3d_torus_initial_state
   )

set(3d_torus_driver_FLAGS COMPONENT_DEPENDENCIES octopus)
set(3d_torus_initial_state_FLAGS COMPONENT_DEPENDENCIES octopus)

if(VISIT_DRIVER_FOUND AND VISIT_LIBSIM_FOUND)
  set(applications
      ${applications}
      3d_torus_realtime_viz
     )
  set(3d_torus_realtime_viz_FLAGS
      COMPONENT_DEPENDENCIES octopus visit_simulation)
endif()

foreach(application ${applications})
  set(sources ${application}.cpp)

  source_group("Source Files" FILES ${sources})

  add_hpx_executable(${application} SOURCES ${sources} ${${application}_FLAGS}
                     FOLDER "Applications")

  add_hpx_pseudo_target(applications.${application})
  add_hpx_pseudo_dependencies(applications applications.${application})
  add_hpx_pseudo_dependencies(applications.${application} ${application}_exe)
endforeach()

