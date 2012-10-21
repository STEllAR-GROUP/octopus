////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(E8D1FEAF_F0FB_47CA_81E6_6EC657FF392F)
#define E8D1FEAF_F0FB_47CA_81E6_6EC657FF392F

#include <octopus/assert.hpp>

#include <ostream>

#include <boost/serialization/split_free.hpp>

namespace octopus
{

enum axis
{
    invalid_axis = -1
  , x_axis       = 0
  , y_axis       = 1
  , z_axis       = 2
};

inline std::ostream& operator<<(std::ostream& os, axis f)
{
    switch (f)
    {
        case x_axis: os << "x_axis"; break; 
        case y_axis: os << "y_axis"; break; 
        case z_axis: os << "z_axis"; break; 
        default: os << "invalid_axis"; break; 
    }
    return os;
}

}

///////////////////////////////////////////////////////////////////////////////
namespace boost { namespace serialization
{
    template <typename Archive>
    void save(Archive& ar, octopus::axis const& k, const unsigned int)
    {
        boost::uint8_t tmp(k);
        ar & tmp; 
    }

    template <typename Archive>
    void load(Archive& ar, octopus::axis& k, const unsigned int)
    {
        boost::uint8_t tmp;
        ar & tmp; 
        k = tmp; 
    }
}}

BOOST_SERIALIZATION_SPLIT_FREE(octopus::axis);

#endif // OCTOPUS_E8D1FEAF_F0FB_47CA_81E6_6EC657FF392F

