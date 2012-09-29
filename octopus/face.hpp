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

namespace octopus
{

enum face
{
    XL = 0, ///< X-lower
    XU = 1, ///< X-upper
    YL = 2, ///< Y-lower
    YU = 3, ///< Y-upper
    ZL = 4, ///< Z-lower
    ZU = 5  ///< Z-upper
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
    }
}

inline std::ostream& operator<<(std::ostream& os, face f)
{
    switch (f)
    {
        case XL: os << "XL";
        case XU: os << "XU";
        case YL: os << "YL";
        case YU: os << "YU";
        case ZL: os << "ZL";
        case ZU: os << "ZU"; 
    }
    return os;
}

}

#endif // OCTOPUS_DE82DE29_3AEC_4D0E_A2FD_AF7C424F9080

