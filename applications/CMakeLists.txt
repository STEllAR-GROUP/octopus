# Copyright (c) 2011-2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(subdirs
    3d_torus
    sod_shock_tube
    rayleigh_taylor
   )

foreach(subdir ${subdirs})
  add_subdirectory(${subdir})
endforeach()

