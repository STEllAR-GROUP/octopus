////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_266BB67A_CF59_461F_B534_6FA6ACBB59F7)
#define OCTOPUS_266BB67A_CF59_461F_B534_6FA6ACBB59F7

#include <octopus/config.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/version.hpp>

#include <iostream>

#define OCTOPUS_CONFIG_DATA_VERSION 0x01

// TODO: This is specific to the euler code, make it more general after SC.

namespace octopus
{

// NOTE (to self): When you add configuration parameters, don't forget to add
// them to configuration_from_init().
// NOTE: This is kept an aggregate for simplicity. Default values are set by
// configuration_from_init().
struct config_data
{
    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Dominic's original code.
    //
    // Format:
    //     ///< description
    //     type name; // original name, original value, notes

    ///< The dimensional size of the grid including ghost zones. I need
    ///  clarification on whether this is the size of the entire grid, or the
    ///  number of points that each grid node represents. 
    boost::uint64_t dimensional_size; // GNX, 8+2*bw, TODO: validate min/max 

    ///////////////////////////////////////////////////////////////////////////
    // "My" parameters (stuff not in Dominic's code).

    ///< This is the number of timesteps in the future that the code should try
    ///  to predict. This is also the number of timesteps alive at any given
    ///  time after ramp up of the code. 
    boost::uint64_t temporal_prediction_limit; 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & dimensional_size;
    }
};

OCTOPUS_EXPORT std::ostream& operator<<(
    std::ostream& os
  , config_data const& cfg
    );

OCTOPUS_EXPORT config_data config_from_ini();

}

BOOST_CLASS_VERSION(octopus::config_data, OCTOPUS_CONFIG_DATA_VERSION)
BOOST_CLASS_TRACKING(octopus::config_data, boost::serialization::track_never)

#endif // OCTOPUS_266BB67A_CF59_461F_B534_6FA6ACBB59F7

