////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E)
#define OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E

#include <boost/cstdint.hpp>

// TODO: Stick one of these into engine_server.

namespace octopus
{

/// A table that does... science. This is a vtable defining the physics through
/// a set constants and serializable functions. Note that the constants in this
/// class are distinctly different from parameters in config_data. Constants
/// in a science table represent values that vary from problem to problem, but
/// are fixed within the scope of each problem. 
// NOTE: Aggregate for laziness.
struct science_table
{
    boost::uint64_t ghost_zone_width; /// The width of ghost zones on all sides,
                                      /// measured in number of grid points. 
    boost::uint64_t state_size;  /// Number of doubles needed for state for
                                 /// each discrete value on the grid.
 
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & ghost_zone_width;
        ar & state_size;
    }
};

}

#endif // OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E

