////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_A24A238E_0CA5_4002_8A9A_2923443578DD)
#define OCTOPUS_A24A238E_0CA5_4002_8A9A_2923443578DD

#include <boost/fusion/include/define_struct.hpp>

#include <hpx/util/serialize_sequence.hpp>

BOOST_FUSION_DEFINE_STRUCT(
    (octopus), timestep_prediction,
    (double, next_dt)
    (double, future_dt))

namespace boost { namespace serialization
{

template <typename Archive>
void serialize(Archive& ar, octopus::timestep_prediction& p, unsigned int)
{
    hpx::util::serialize_sequence(ar, p);
}

}}

#endif // OCTOPUS_A24A238E_0CA5_4002_8A9A_2923443578DD

