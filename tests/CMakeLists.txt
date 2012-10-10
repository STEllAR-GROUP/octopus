# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(tests
    child_index
    octree_create_child
    vector3d
   )

foreach(test ${tests})
  set(sources ${test}.cpp)

  source_group("Source Files" FILES ${sources})

  add_hpx_executable(${test}_test SOURCES ${sources} ${${test}_FLAGS}
                     DEPENDENCIES octopus_component
                     FOLDER "Tests")

  add_hpx_pseudo_target(tests.${test})
  add_hpx_pseudo_dependencies(tests tests.${test})
  add_hpx_pseudo_dependencies(tests.${test} ${test}_test_exe)
endforeach()

