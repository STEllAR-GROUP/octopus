////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_BAD2F15A_83F9_4746_A37A_EF33AE76B481)
#define OCTOPUS_BAD2F15A_83F9_4746_A37A_EF33AE76B481

#include <octopus/config.hpp>

#include <ostream>

namespace octopus
{

enum module_role
{
    invalid_role     = -1,
    node_distributor = 0,
    last_role
};

inline std::ostream& operator<<(std::ostream& os, module_role r)
{
    switch (r)
    {
        case node_distributor: os << "node_distributor"; break; 
        default: os << "invalid_role"; break; 
    }
    return os;
}

inline std::string default_module(module_role r)
{
    switch (r)
    {
        case node_distributor: return "single_locality_distributor"; 
        default: return "invalid_role"; 
    }
}

struct OCTOPUS_EXPORT module_base
{
    virtual ~module_base() {}

    virtual module_role role() const = 0;
};

}

#endif // OCTOPUS_BAD2F15A_83F9_4746_A37A_EF33AE76B481

