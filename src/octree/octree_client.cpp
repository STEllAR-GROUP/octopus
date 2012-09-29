////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/lcos/future.hpp>
#include <hpx/async.hpp>

#include <octopus/octree/octree_server.hpp>

namespace octopus
{

void octree_client::create(hpx::id_type const& locality)
{
    OCTOPUS_ASSERT_MSG(locality.get_msb() & 0xFF
                     , "target GID is not a locality");

    const hpx::components::component_type t
        = hpx::components::get_component_type<octopus::octree_server>();

    hpx::components::runtime_support rts(locality);
    gid_ = rts.create_component(t, 1);
}

void octree_client::create_child(child_index idx)
{
    create_child_async(idx).get(); 
}

hpx::future<void> octree_client::create_child_async(child_index idx)
{
    return hpx::async<octree_server::create_child_action>(gid_, idx);
}

}

