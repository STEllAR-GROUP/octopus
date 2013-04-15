////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_17914F5A_C09A_42A4_9819_0641B54EBF15)
#define OCTOPUS_17914F5A_C09A_42A4_9819_0641B54EBF15

#include <octopus/config.hpp>
#include <octopus/science/science_table.hpp>

#include <boost/program_options.hpp>

// You must define at least one of these two functions. You can define both if
// you wish.

/// octopus_define_problem is called after the AMR framework has been
/// initialized. It is passed an octopus::science_table filled with default
/// values. 
extern "C" 
void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    ) OCTOPUS_WEAK;

/// octopus_main is called after the AMR framework has been initialized and the
/// problem defined.
extern "C"
int octopus_main(
    boost::program_options::variables_map& vm
    ) OCTOPUS_WEAK;

#endif // OCTOPUS_17914F5A_C09A_42A4_9819_0641B54EBF15

