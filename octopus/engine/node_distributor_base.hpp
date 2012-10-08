////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_AD19F7EF_D49E_4A45_A083_12FA5ED535E7)
#define OCTOPUS_AD19F7EF_D49E_4A45_A083_12FA5ED535E7

#include <hpx/lcos/future.hpp>

#include <octopus/array1d.hpp>
#include <octopus/engine/module.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT node_distributor_base : module_base
{
    virtual ~node_distributor_base() {}

    virtual hpx::future<hpx::id_type, hpx::naming::gid_type> create_async(
        boost::uint64_t level
      , array1d<boost::uint64_t, 3> const& location
      , hpx::id_type const& hint = hpx::naming::invalid_id
        ) = 0;

    virtual hpx::id_type create(
        boost::uint64_t level
      , array1d<boost::uint64_t, 3> const& location
      , hpx::id_type const& hint = hpx::naming::invalid_id
        )
    {
        return create_async(level, location, hint).get();  
    }

    module_role role() const
    {
        return node_distributor; 
    }
};

}

#endif // OCTOPUS_AD19F7EF_D49E_4A45_A083_12FA5ED535E7

