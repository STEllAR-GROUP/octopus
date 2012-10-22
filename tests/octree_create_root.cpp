////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_fwd.hpp>
#include <hpx/util/lightweight_test.hpp>

#include <octopus/science.hpp>
#include <octopus/driver.hpp>
#include <octopus/octree/octree_client.hpp>
#include <octopus/octree/checkout_state.hpp>

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    root.create_root(hpx::find_here(), octopus::octree_init_data());

    {
        boost::array<octopus::octree_client, 6> siblings = root.get_siblings();
    
        for ( octopus::face i = octopus::XL
            ; i < octopus::invalid_face
            ; i = octopus::face(boost::uint8_t(i + 1)))
        {
            HPX_TEST(!siblings[i].real());
            HPX_TEST_EQ(siblings[i].kind(), octopus::physical_boundary);
        }
    }

    return hpx::util::report_errors();
}

