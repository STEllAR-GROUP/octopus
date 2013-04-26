////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/operators/boost_array_arithmetic.hpp>
#include <octopus/science/reconstruction.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/math.hpp>

namespace octopus
{

void minmod_reconstruction::operator()(
    std::vector<state> const& q0
  , std::vector<state>& ql
  , std::vector<state>& qr
    ) const
{
    boost::uint64_t const gnx = config().grid_node_length;

    std::vector<state> slope(gnx, state());

    for (boost::uint64_t i = 1; i < gnx - 1; ++i)
    {
        // up = q0[i + 1] - q0[i]
        state up  = q0[i + 1];
                            up -= q0[i];

        // um = q0[i] - q0[i - 1]
        state um  = q0[i];
                            um -= q0[i - 1];

        slope[i] = minmod_theta(up, um, theta_);
    }

    for (boost::uint64_t i = 2; i < gnx - 1; ++i)
    {
        // ql[i] = q0[i - 1] + slope[i - 1] / 2.0
        ql[i]  = q0[i - 1];
        ql[i] += slope[i - 1] / 2.0;

        //qr[i] = q0[i] - slope[i] / 2.0;
        qr[i]  = q0[i];
        qr[i] -= slope[i] / 2.0;
    }
}

}

