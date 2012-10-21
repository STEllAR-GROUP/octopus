////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E)
#define OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E

#include <octopus/face.hpp>

#include <boost/cstdint.hpp>

// NOTE: (to self) Don't forgot to update default_science_table when
// science_table is updated.
// FIXME: This should probably live in another header.

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
    boost::uint64_t state_size;  /// Number of doubles needed for state for
                                 /// each discrete value (zone) on the grid.

    ///////////////////////////////////////////////////////////////////////////
    // State getters. FIXME: Abstract these away entirely into application code.
    // Any code in the AMR driver that needs to be aware of what's in the state
    // should be in the science table instead. 
    hpx::util::function<
        double const&(std::vector<double> const&)
    > rho; 

    /// Defines the physical boundaries. Returns true if a face at a location in
    /// the octree and a particular level of refinement is a physical boundary. 
    hpx::util::function<
        bool(
            boost::array<boost::int64_t, 3> const&
          , face
          , boost::uint64_t
            )
    > physical_boundaries; 

    /// The reconstruction scheme.
    // FIXME: How did an std::vector<std::vector<> > sneak into this code? Use
    // vector3d instead.
    hpx::util::function<
        void(
            std::vector<std::vector<double> > const&
          , std::vector<std::vector<double> >&
          , std::vector<std::vector<double> >&
            )
    > reconstruction; 

    boost::uint64_t ghost_zone_width; /// The width of ghost zones on all sides,
                                      /// measured in number of grid points. 

    hpx::util::function<
        void(
            boost::array<double, 3> const&
          , face
            )
    > enforce_outflow; 

    hpx::util::function<
        void(std::vector<double>&)
    > reflect_on_x; 

    hpx::util::function<
        void(std::vector<double>&)
    > reflect_on_y; 

    hpx::util::function<
        void(std::vector<double>&)
    > reflect_on_z; 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & state_size;
        ar & rho; 

        ar & physical_boundaries;
        ar & reconstruction;
        ar & ghost_zone_width;
        ar & enforce_outflow;

        ar & reflect_on_x;
        ar & reflect_on_y;
        ar & reflect_on_z;
    }
};

OCTOPUS_EXPORT science_table default_science_table();

}

#endif // OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E

