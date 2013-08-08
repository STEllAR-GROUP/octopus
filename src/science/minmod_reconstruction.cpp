////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/science/minmod_reconstruction.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/math.hpp>

namespace octopus
{

// FIXME: Put this in a header?
void minmod_reconstruction::operator()(
/*
    vector2d<double> const& q0
  , vector2d<double>& ql
  , vector2d<double>& qr
*/
    std::vector<state> const& q0
  , std::vector<state>& ql
  , std::vector<state>& qr
    ) const
{
    boost::uint64_t const gnx = config().grid_node_length;

//    vector2d<double> slope(gnx);
    std::vector<state> slope(gnx);

    for (boost::uint64_t i = 1; i < gnx - 1; ++i)
    {
        state up = q0[i + 1] - q0[i];
        state um = q0[i] - q0[i - 1];
        slope[i] = minmod_theta(up, um, theta_);
    }

    for (boost::uint64_t i = 2; i < gnx - 1; ++i)
    {
        ql[i] = q0[i - 1] + (slope[i - 1] / 2.0);
        qr[i] = q0[i] - (slope[i] / 2.0);
    }
}

}

