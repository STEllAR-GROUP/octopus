////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_4A14DCB9_DB45_4F44_BFD2_4DFE2E53F6D6)
#define OCTOPUS_4A14DCB9_DB45_4F44_BFD2_4DFE2E53F6D6

#include <octopus/assert.hpp>
#include <octopus/face.hpp>
#include <octopus/array.hpp>
#include <octopus/trivial_serialization.hpp>

namespace octopus
{

struct physical_boundaries_at_zero : trivial_serialization
{
    bool operator()(
        array<boost::int64_t, 3> const& location
      , face f
      , boost::uint64_t level
        ) const
    {
        switch (f)
        {
            ///////////////////////////////////////////////////////////////////
            // X-axis.
            case XL:
                return location[0] == 0;
            case XU:
                return location[0] == ((1 << level) - 1);
    
            ///////////////////////////////////////////////////////////////////
            // Y-axis.
            case YL:
                return location[1] == 0;
            case YU:
                return location[1] == ((1 << level) - 1);
    
            ///////////////////////////////////////////////////////////////////
            // Z-axis.
            case ZL:
                return location[2] == 0;
            case ZU:
                return location[2] == ((1 << level) - 1);

            default:
                break;
        } 
    
        OCTOPUS_ASSERT(false);
        return false;
    }
};

}

#endif // OCTOPUS_4A14DCB9_DB45_4F44_BFD2_4DFE2E53F6D6

