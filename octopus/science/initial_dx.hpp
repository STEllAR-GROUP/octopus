////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_7AD74C5D_0848_4DC1_98A4_8687A5D6CF00)
#define OCTOPUS_7AD74C5D_0848_4DC1_98A4_8687A5D6CF00

#include <octopus/engine/engine_interface.hpp>
#include <octopus/trivial_serialization.hpp>

namespace octopus
{

struct initial_dx : trivial_serialization
{
    double operator()() const
    {
        double const grid_dim = octopus::config().spatial_domain;
        boost::uint64_t const gnx = octopus::config().grid_node_length;
        boost::uint64_t const bw = octopus::science().ghost_zone_length; 

        return (2.0 * grid_dim / double(gnx - 2 * bw));
    }
};

}

#endif // OCTOPUS_7AD74C5D_0848_4DC1_98A4_8687A5D6CF00

