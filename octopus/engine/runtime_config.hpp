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
// them to config_from_init().
// NOTE: This is kept an aggregate for simplicity. Default values are set by
// config_from_init().
// NOTE: Users should never copy this.
struct config_data
{
    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Dominic's original code.
    //
    // Format:
    //     ///< description
    //     type name; // original name, original value, notes

    boost::uint64_t max_refinement_level; 

    ///< The edge length of each grid node (units == discrete points aka zones) 
    ///  including ghost zones.
    /// NOTE: Needs confirmation. 
    // NOTE: This MUST be a power of 2.
    boost::uint64_t grid_node_length; // GNX, 8+2*BW, TODO: validate min/max 

    ///< The spatial edge length of the problem NOTE: Confirm.
    double spatial_domain; // GRID_DIM, 1.5e-4 (Zach) and 1.0 (Dominic) 

    ///< The "physics" distance between zones in the root node/most coarse
    ///  refinement level (aka level 0). 
    /// NOTE: Needs confirmation. 
    /// NOTE: Rename to initial_space_step 
    double initial_spatial_step; // h0, (2*GRID_DIM/double(GNX-2*BW)),
                                 // TODO: validate min/max

    ///< Order of (TVD) Runge Kutta used to solve the PDE.
    boost::uint16_t runge_kutta_order;

    ///< Reflection control for each axis. TODO: error handling if no reflection
    ///  function is available for a particular axis.
    bool x_reflect;
    bool y_reflect;
    bool z_reflect;
 
    ///////////////////////////////////////////////////////////////////////////
    // "My" parameters (stuff not in Dominic's code).

    ///< This is the gap, in timesteps, between a timestep, and the timestep
    ///  that computed its timestep size. For example, if this is 10, then,
    ///  timestep 1 computes the timestep size for timestep 11, and timestep
    ///  2 computes the the timestep size for timestep 12, etc. This is also the
    ///  number of timesteps alive at any given time after ramp up of the code. 
    boost::uint64_t temporal_prediction_gap; 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & max_refinement_level;
        ar & grid_node_length;
        ar & spatial_domain;
        ar & initial_spatial_step;
        ar & runge_kutta_order;
        ar & x_reflect;
        ar & y_reflect;
        ar & z_reflect;
        ar & temporal_prediction_gap;
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

