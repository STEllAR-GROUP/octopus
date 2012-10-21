////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/operators/std_vector_arithmetic.hpp>
#include <octopus/science/reconstruction.hpp>
#include <octopus/engine/engine_interface.hpp>

namespace octopus
{

void minmod_reconstruction::operator()(
    std::vector<std::vector<double> > const& q0
  , std::vector<std::vector<double> >& ql
  , std::vector<std::vector<double> >& qr
    ) const
{
    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const gnx = config().grid_node_length;

    std::vector<std::vector<double> > slope(gnx, std::vector<double>(ss));

    using namespace octopus::operators;

    for (boost::uint64_t i = 1; i < gnx - 1; ++i)
    {
        // up = q0[i + 1] - q0[i]
        std::vector<double> up  = q0[i + 1];
                            up -= q0[i];

        // um = q0[i] - q0[i - 1]
        std::vector<double> um  = q0[i];
                            um -= q0[i - 1];

        for (boost::uint64_t l = 0; l < ss; ++l)
        {
            if (up[l] * um[l] > 0.0)
                slope[i][l] = 2.0 * up[l] * um[l] / (up[l] + um[l]);
            else
                slope[i][l] = 0.0;
        }
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

