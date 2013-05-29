# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(tests
    global_variable
   )

foreach(application ${tests})
  set(sources ${application}.cpp)

  source_group("Source Files" FILES ${sources})

  add_hpx_executable(${application} SOURCES ${sources} ${${application}_FLAGS}
    FOLDER "Tests")

  add_hpx_pseudo_target(tests.${application})
  add_hpx_pseudo_dependencies(tests tests.${application})
  add_hpx_pseudo_dependencies(tests.${application} ${application}_exe)
endforeach()

