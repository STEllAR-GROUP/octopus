////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_E73B96D9_CF42_4331_A4D5_2BE254FD289C)
#define OCTOPUS_E73B96D9_CF42_4331_A4D5_2BE254FD289C

#include <octopus/array1d.hpp>
#include <octopus/vector3d.hpp>

#include <boost/cstdint.hpp>

namespace octopus
{

// FIXME: Tactic needed for constructing the init data for root.
// NOTE: Aggregate for laziness.
struct octree_init_data
{
    boost::uint64_t                 level;
    array1d<boost::int64_t, 3>      location;
    double                          dx;
    double                          time; 
    array1d<boost::int64_t, 3>      offset; 
    array1d<double, 3>              origin;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & level;
        ar & location;
        ar & dx;
        ar & time;
        ar & offset;
        ar & origin;
    }
};

}

#endif // OCTOPUS_E73B96D9_CF42_4331_A4D5_2BE254FD289C

