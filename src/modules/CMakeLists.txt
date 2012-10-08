# Copyright (c) 2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(modules
    single_locality_distributor
    #round_robin_distributor
   )

foreach(module ${modules})
  add_hpx_component(${module}
    SOURCES ${module}.cpp
    FOLDER "Modules"
    ESSENTIAL
    OUTPUT_SUFFIX lib/hpx/octopus
    INSTALL_SUFFIX lib/hpx/octopus)
endforeach()


