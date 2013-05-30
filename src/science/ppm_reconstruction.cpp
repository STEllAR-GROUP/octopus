////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/science/ppm_reconstruction.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/math.hpp>

namespace octopus
{

// FIXME: Put this in a header?
void ppm_reconstruction::operator()(
    std::vector<state> const& q0
  , std::vector<state>& ql
  , std::vector<state>& qr
    ) const
{
    boost::uint64_t const gnx = config().grid_node_length;

    std::vector<state> slope(gnx, state());

    for (boost::uint64_t i = 1; i < gnx - 1; ++i)
    {
        state up = q0[i + 1] - q0[i];
        state um = q0[i] - q0[i - 1];
        slope[i] = minmod_theta(up, um, 2.0);
    }

    for (boost::uint64_t i = 1; i < gnx - 2; ++i)
    {
        ql[i] = (q0[i] + q0[i + 1]) * 0.5;
        for (boost::uint64_t l = 0; l < slope[i].size(); ++l)
            ql[i][l] += (slope[i][l] - slope[i + 1][l]) * (1.0 / 6.0);
        qr[i] = ql[i - 1];
    }

    for (boost::uint64_t i = 2; i < gnx - 2; ++i)
    {
        for (boost::uint64_t l = 0; l < slope[i].size(); ++l)
        {
            double const t0 = ql[i][l] - qr[i][l];
            double const t1 = ql[i][l] + qr[i][l];

            if ((ql[i][l] - q0[i][l]) * (q0[i][l] - qr[i][l]) <= 0.0)
                ql[i][l] = qr[i][l] = q0[i][l];
            else if (t0 * (q0[i][l] - 0.5 * t1) > (1.0 / 6.0) * t0 * t0)
                qr[i][l] = 3.0 * q0[i][l] - 2.0 * ql[i][l];
            else if (-(1.0 / 6.0) * t0 * t0 > t0 * (q0[i][l] - 0.5 * t1))
                ql[i][l] = 3.0 * q0[i][l] - 2.0 * qr[i][l];
        }
    }

    for (boost::uint64_t i = gnx - 3; i > 2; --i)
        ql[i] = ql[i - 1];
}

}

