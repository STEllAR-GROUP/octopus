////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/lcos/future.hpp>
#include <hpx/async.hpp>
#include <hpx/apply.hpp>
#include <hpx/runtime/components/runtime_support.hpp>

#include <octopus/octree/octree_server.hpp>

namespace octopus
{

void octree_client::create(
    hpx::id_type const& locality
    )
{
    OCTOPUS_ASSERT_FMT_MSG(locality.get_msb() & 0xFF,
                           "target is not a locality, gid(%1%)",
                           locality);

    const hpx::components::component_type t
        = hpx::components::get_component_type<octopus::octree_server>();

    hpx::components::runtime_support rts(locality);
    gid_ = rts.create_component(t, 1);
}

hpx::future<hpx::id_type, hpx::naming::gid_type> octree_client::create_async(
    hpx::id_type const& locality
    ) const
{
    OCTOPUS_ASSERT_FMT_MSG(locality.get_msb() & 0xFF,
                           "target is not a locality, gid(%1%)",
                           locality);

    const hpx::components::component_type t
        = hpx::components::get_component_type<octopus::octree_server>();

    hpx::components::runtime_support rts(locality);
    return rts.create_component_async(t, 1);
}

void octree_client::create_child(
    child_index kid
    )
{
    create_child_async(kid).get(); 
}

hpx::future<void> octree_client::create_child_async(
    child_index kid
    )
{
    return hpx::async<octree_server::create_child_action>(gid_, kid);
}

void octree_client::set_sibling(
    boost::uint8_t f
  , octree_client const& sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::async<octree_server::set_sibling_action>(gid_, f, sib).get();
}

void octree_client::tie_sibling(
    boost::uint8_t target_f
  , octree_client const& target_sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::async<octree_server::tie_sibling_action>
        (gid_, target_f, target_sib).get();
}

void octree_client::tie_sibling_push(
    boost::uint8_t target_f
  , octree_client const& target_sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::apply<octree_server::tie_sibling_action>
        (gid_, target_f, target_sib);
}

void octree_client::set_sibling_push(
    boost::uint8_t f
  , octree_client const& sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::apply<octree_server::set_sibling_action>(gid_, f, sib);
}

void octree_client::set_child_sibling(
    child_index kid
  , boost::uint8_t f
  , octree_client const& sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::async<octree_server::set_child_sibling_action>
        (gid_, kid, f, sib).get();
}

void octree_client::set_child_sibling_push(
    child_index kid
  , boost::uint8_t f
  , octree_client const& sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::apply<octree_server::set_child_sibling_action>(gid_, kid, f, sib);
}

void octree_client::tie_child_sibling(
    child_index target_kid
  , boost::uint8_t target_f
  , octree_client const& target_sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::async<octree_server::set_child_sibling_action>
        (gid_, target_kid, target_f, target_sib).get();
}

void octree_client::tie_child_sibling_push(
    child_index target_kid
  , boost::uint8_t target_f
  , octree_client const& target_sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::apply<octree_server::set_child_sibling_action>
        (gid_, target_kid, target_f, target_sib);
}

}

