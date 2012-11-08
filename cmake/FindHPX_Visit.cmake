# Copyright (c) 2011-2012 Bryce Adelstein-Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

if(NOT HPX_FINDPACKAGE_LOADED)
  include(HPX_FindPackage)
endif()

if(NOT HPX_FINDPROGRAM_LOADED)
  include(HPX_FindProgram)
endif()

hpx_find_program(VISIT_DRIVER
  PROGRAMS visit
  PROGRAM_PATHS bin)

hpx_find_package(VISIT_LIBSIM
  LIBRARIES simV2 
  LIBRARY_PATHS lib 
  HEADERS VisItControlInterface_V2.h 
  HEADER_PATHS include)

