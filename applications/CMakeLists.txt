# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(applications
    3d_torus
   )

foreach(test ${applications})
  set(sources ${test}.cpp)

  source_group("Source Files" FILES ${sources})

  add_hpx_executable(${test}_test SOURCES ${sources} ${${test}_FLAGS}
                     DEPENDENCIES octopus_component
                     FOLDER "Applications")

  add_hpx_pseudo_target(applications.${test})
  add_hpx_pseudo_dependencies(applications applications.${test})
  add_hpx_pseudo_dependencies(applications.${test} ${test}_test_exe)
endforeach()

