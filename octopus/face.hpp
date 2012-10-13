////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_DE82DE29_3AEC_4D0E_A2FD_AF7C424F9080)
#define OCTOPUS_DE82DE29_3AEC_4D0E_A2FD_AF7C424F9080

#include <ostream>

#include <boost/serialization/split_free.hpp>

namespace octopus
{

enum face
{
                      // direction
    XL = 0, // X-lower: (-1,  0,  0) 
    XU = 1, // X-upper: (+1,  0,  0) 
    YL = 2, // Y-lower: ( 0, -1,  0)
    YU = 3, // Y-upper: ( 0, +1,  0) 
    ZL = 4, // Z-lower: ( 0,  0, -1)
    ZU = 5, // Z-upper: ( 0,  0, +1)

    out_of_bounds
};

inline face invert(face f)
{
    switch (f)
    {
        case XL: return XU;
        case XU: return XL;
        case YL: return YU;
        case YU: return YL;
        case ZL: return ZU;
        case ZU: return ZL; 
        default: break;
    }

    OCTOPUS_ASSERT_MSG(false, "attempt to invert out-of-bounds face"); 
    return out_of_bounds;
}

inline std::ostream& operator<<(std::ostream& os, face f)
{
    switch (f)
    {
        case XL: os << "XL"; break; 
        case XU: os << "XU"; break; 
        case YL: os << "YL"; break; 
        case YU: os << "YU"; break; 
        case ZL: os << "ZL"; break;  
        case ZU: os << "ZU"; break; 
        default: os << "out_of_bounds"; break; 
    }
    return os;
}

}

///////////////////////////////////////////////////////////////////////////////
namespace boost { namespace serialization
{
    template <typename Archive>
    void save(Archive& ar, octopus::face const& k, const unsigned int)
    {
        boost::uint8_t tmp(k);
        ar & tmp; 
    }

    template <typename Archive>
    void load(Archive& ar, octopus::face& k, const unsigned int)
    {
        boost::uint8_t tmp;
        ar & tmp; 
        k = tmp; 
    }
}}

BOOST_SERIALIZATION_SPLIT_FREE(octopus::face);

#endif // OCTOPUS_DE82DE29_3AEC_4D0E_A2FD_AF7C424F9080

