////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/driver.hpp>
#include <octopus/octree/octree_server.hpp>
#include <octopus/engine/engine_interface.hpp>

#include <boost/program_options.hpp>

void octopus_define_problem(
    boost::program_options::variables_map& 
  , octopus::science_table&
    )
{
    // No-op.
}

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    octopus::octree_init_data root_data;
    root_data.dx = octopus::science().initial_spacestep();
    root.create_root(hpx::find_here(), root_data);

    root.apply(octopus::science().initialize);

    root.refine();

/* IMPLEMENT
    if (config().enable_output)
        root.output_initial();
    root.walk();
*/

    return 0;
}

