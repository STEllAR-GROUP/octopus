////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/runtime/components/runtime_support.hpp>

#include <typeinfo>

#include <boost/plugin.hpp>

#include <octopus/engine/node_distributor_base.hpp>
#include <octopus/octree/octree_server.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT single_locality_distributor : node_distributor_base
{
    hpx::future<hpx::id_type, hpx::naming::gid_type> create_async(   
        boost::uint64_t level
      , array1d<boost::uint64_t, 3> const& location
      , hpx::id_type const& hint = hpx::naming::invalid_id
        )
    {
        hpx::components::runtime_support rts(hpx::find_here());
        return rts.create_component_async<octopus::octree_server>
            (level, location);
    }
};

}

BOOST_PLUGIN_EXPORT(octopus,
                    octopus::module_base,
                    octopus::single_locality_distributor,
                    single_locality_distributor,
                    module);

BOOST_PLUGIN_EXPORT_LIST(octopus, module);

