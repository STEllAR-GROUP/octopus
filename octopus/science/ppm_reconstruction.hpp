////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_125247F0_F7B5_42E1_B3DD_957000B62E58)
#define OCTOPUS_125247F0_F7B5_42E1_B3DD_957000B62E58

#include <octopus/config.hpp>
#include <octopus/state.hpp>
#include <octopus/trivial_serialization.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT ppm_reconstruction : trivial_serialization
{
    enum { ghost_zone_width = 3 };

    void operator()(
        std::vector<state> const& q0
      , std::vector<state>& ql
      , std::vector<state>& qr
        ) const;
};

}

#endif // OCTOPUS_125247F0_F7B5_42E1_B3DD_957000B62E58

