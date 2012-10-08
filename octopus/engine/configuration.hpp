////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_266BB67A_CF59_461F_B534_6FA6ACBB59F7)
#define OCTOPUS_266BB67A_CF59_461F_B534_6FA6ACBB59F7

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/version.hpp>

#define OCTOPUS_CONFIGURATION_VERSION 0x01

namespace octopus
{

struct configuration
{
                            // In the original code:
    boost::uint64_t bw;     // 2
    boost::uint64_t gnx;    // 8+2*bw 
    boost::uint64_t nzones; // (gnx-2*bw)*(gnx-2*bw)*(gnx-2*bw) 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & bw;
        ar & gnx;
        ar & nzones;
    }
};

}

BOOST_CLASS_VERSION(octopus::configuration, OCTOPUS_CONFIGURATION_VERSION)
BOOST_CLASS_TRACKING(octopus::configuration, boost::serialization::track_never)

#endif // OCTOPUS_266BB67A_CF59_461F_B534_6FA6ACBB59F7

