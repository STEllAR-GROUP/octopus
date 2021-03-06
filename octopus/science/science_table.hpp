////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E)
#define OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E

#include <octopus/octree/octree_server.hpp>
#include <octopus/science/refinement_criteria.hpp>
#include <octopus/science/dt_prediction.hpp>
#include <octopus/vector2d.hpp>
#include <octopus/io/writer.hpp>
#include <octopus/face.hpp>

#include <boost/cstdint.hpp>

// NOTE: (to self) Don't forgot to update default_science_table when
// science_table is updated.

#define OCTOPUS_SCIENCE_TABLE_VERSION 0x01

namespace octopus
{

/// A table that does... science. This is a vtable defining the physics through
/// a set constants and serializable functions. Note that the constants in this
/// class are distinctly different from parameters in config_data. Constants
/// in a science table represent values that vary from problem to problem, but
/// are fixed within the scope of each problem. 
// NOTE: Aggregate for laziness.
// NOTE: Users should never copy this.
struct science_table
{
    /// Defines the physical boundaries. Returns true if a face at a location in
    /// the octree and a particular level of refinement is a physical boundary. 
    hpx::util::function<
        bool(
            array<boost::int64_t, 3> const&
          , face
          , boost::uint64_t
            )
    > physical_boundaries; 

    /// The reconstruction scheme.
    hpx::util::function<
        void(
/*
            vector2d<double> const&
          , vector2d<double>&
          , vector2d<double>&
*/
            std::vector<state> const&
          , std::vector<state>&
          , std::vector<state>&
            )
    > reconstruct; 

    boost::uint64_t ghost_zone_length; /// The width of ghost zones on all sides,
                                      /// measured in number of grid points. 

    hpx::util::function<
        void(octree_server&)
    > initialize;

    hpx::util::function<
        void(
            octree_server&
          , state& 
          , array<double, 3> const& 
          , face 
        )
    > enforce_outflow; 

    hpx::util::function<
        void(
            state& 
          , array<double, 3> const&
            )
    > enforce_limits; 

    hpx::util::function<
        void(state&)
    > reflect_z; 

    hpx::util::function<
        double(
            octree_server&
          , state const& 
          , array<double, 3> const& 
          , axis
            )
    > max_eigenvalue; 

    ///< The distance between zones in the root node/most coarse / refinement
    ///  level (aka level 0). 
    // h0, (2*GRID_DIM/double(GNX-2*BW)), TODO: validate min/max
    hpx::util::function<double()> initial_dx;

    hpx::util::function<
        double(
            octree_server& ///< Root
            )
    > initial_dt;

    hpx::util::function<
        /// (timestep N + 1 size, timestep N + gap size)
        dt_prediction(
            octree_server& ///< Root
            )
    > predict_dt;

    hpx::util::function<
        void(
            state& 
          , array<double, 3> const&
            )
    > conserved_to_primitive; 

    hpx::util::function<
        void(
            state& 
          , array<double, 3> const&
            )
    > primitive_to_conserved; 

    hpx::util::function<
        state(
            octree_server&
          , state const& 
          , array<double, 3> const&
            )
    > source; 

    hpx::util::function<
        state(
            octree_server&
          , state&
          , array<double, 3> const& 
          , array<boost::uint64_t, 3> const& 
          , axis
            )
    > flux; 

    hpx::util::function<
        hpx::id_type(
            octree_init_data const& init
          , std::vector<hpx::id_type> const& localities
            )
    > distribute; 

    // FIXME: Rename to refine.
    refinement_criteria refine_policy;

    writer output;

    science_table()
     : physical_boundaries()
     , reconstruct()
     , ghost_zone_length(0)
     , initialize()
     , enforce_outflow()
     , enforce_limits()
     , reflect_z()
     , max_eigenvalue()
     , initial_dx()
     , initial_dt()
     , predict_dt()
     , conserved_to_primitive()
     , primitive_to_conserved()
     , source()
     , flux()
     , distribute()
     , refine_policy()
     , output()
    {}

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & physical_boundaries;
        ar & reconstruct;
        ar & ghost_zone_length;

        ar & initialize;

        ar & enforce_outflow;
        ar & enforce_limits;

        ar & reflect_z;
        ar & max_eigenvalue;

        ar & initial_dx;

        ar & initial_dt;
        ar & predict_dt;

        ar & conserved_to_primitive;
        ar & primitive_to_conserved;

        ar & source;
        ar & flux;

        ar & distribute;
        ar & refine_policy;

        ar & output;
    }
};

OCTOPUS_EXPORT science_table default_science_table();

}

BOOST_CLASS_VERSION(octopus::science_table, OCTOPUS_SCIENCE_TABLE_VERSION)
BOOST_CLASS_TRACKING(octopus::science_table, boost::serialization::track_never)

#endif // OCTOPUS_903E5732_EB7C_4CBC_A13E_24DF0582AF0E

