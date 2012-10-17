////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/child_index.hpp>

namespace octopus
{

child_index invert(face f, child_index const& idx)
{
    child_index idx_inv = idx;

    // Invert 
    switch (f)
    {
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        case XL: // idx_inv = idx + (+1, 0, 0) 
        {
            OCTOPUS_ASSERT(idx.x() == 0);
            idx_inv.set_x(1);
            return idx_inv;
        } 
        case XU: // idx_inv = idx + (-1, 0, 0) 
        {
            OCTOPUS_ASSERT(idx.x() == 1);
            idx_inv.set_x(0);
            return idx_inv;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        case YL: // idx_inv = idx + (0, +1, 0) 
        {
            OCTOPUS_ASSERT(idx.y() == 0);
            idx_inv.set_y(1);
            return idx_inv;
        } 
        case YU: // idx_inv = idx + (0, -1, 0) 
        {
            OCTOPUS_ASSERT(idx.y() == 1);
            idx_inv.set_y(0);
            return idx_inv;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        case ZL: // idx_inv = idx + (0, +1, 0) 
        {
            OCTOPUS_ASSERT(idx.z() == 0);
            idx_inv.set_z(1);
            return idx_inv;
        } 
        case ZU: // idx_inv = idx + (0, -1, 0) 
        {
            OCTOPUS_ASSERT(idx.z() == 1);
            idx_inv.set_z(0);
            return idx_inv;
        }

        default: break;
    }

    OCTOPUS_ASSERT_MSG(false, "face shouldn't be out-of-bounds");
    return child_index();
}

}

