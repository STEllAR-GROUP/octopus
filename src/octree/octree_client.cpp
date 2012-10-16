////////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
void octree_client::create_root(
    hpx::id_type const& locality
  , octree_init_data const& init
    ) 
{
    ensure_real();

    OCTOPUS_ASSERT_FMT_MSG(!(locality.get_msb() & 0xFF),
                           "target is not a locality, gid(%1%)",
                           locality);

    hpx::components::runtime_support rts(locality);
    gid_ = rts.create_component_async<octopus::octree_server>(init).get();
}

void octree_client::create_root(
    hpx::id_type const& locality
  , BOOST_RV_REF(octree_init_data) init
    ) 
{
    ensure_real();

    OCTOPUS_ASSERT_FMT_MSG(!(locality.get_msb() & 0xFF),
                           "target is not a locality, gid(%1%)",
                           locality);

    hpx::components::runtime_support rts(locality);
    gid_ = rts.create_component_async<octopus::octree_server>(init).get();
}

///////////////////////////////////////////////////////////////////////////////
hpx::future<void> octree_client::create_child_async(
    child_index kid
    )
{
    ensure_real();
    return hpx::async<octree_server::create_child_action>(gid_, kid);
}

void octree_client::set_sibling_for_amr_boundary(
    face f
  , octree_client const& sib 
    )
{
    // IMPLEMENT
}

void octree_client::tie_sibling_for_amr_boundary(
    face f
  , octree_client const& sib 
    )
{
    // IMPLEMENT
}

void octree_client::set_sibling_for_physical_boundary(
    face f
  , octree_client const& sib 
    )
{
    // IMPLEMENT
}

void octree_client::tie_sibling_for_physical_boundary(
    face f
  , octree_client const& sib 
    )
{
    // IMPLEMENT
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::set_sibling(
    face f
  , octree_client const& sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));

    if (amr_boundary == kind_)
    {
        set_sibling_for_amr_boundary(f, sib); 
        return;
    }

    else if (physical_boundary == kind_)
    {
        set_sibling_for_physical_boundary(f, sib); 
        return;
    }

    hpx::async<octree_server::set_sibling_action>(gid_, f, sib).get();
}

void octree_client::set_sibling_push(
    face f
  , octree_client const& sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));

    if (amr_boundary == kind_)
    {
        hpx::apply(boost::bind(
            &octree_client::set_sibling_for_amr_boundary, this, _1, _2),
                f, sib);
        return;
    }

    else if (physical_boundary == kind_)
    {
        hpx::apply(boost::bind(
            &octree_client::set_sibling_for_physical_boundary, this, _1, _2),
                f, sib);
        return;
    }

    hpx::apply<octree_server::set_sibling_action>(gid_, f, sib);
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::tie_sibling(
    face target_f
  , octree_client const& target_sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));

    if (amr_boundary == kind_)
    {
        tie_sibling_for_amr_boundary(target_f, target_sib); 
        return;
    }

    else if (physical_boundary == kind_)
    {
        tie_sibling_for_physical_boundary(target_f, target_sib); 
        return;
    }

    hpx::async<octree_server::tie_sibling_action>
        (gid_, target_f, target_sib).get();
}

void octree_client::tie_sibling_push(
    face target_f
  , octree_client const& target_sib
    )
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));

    if (amr_boundary == kind_)
    {
        hpx::apply(boost::bind(
            &octree_client::tie_sibling_for_amr_boundary, this, _1, _2),
                target_f, target_sib);
        return;
    }

    else if (physical_boundary == kind_)
    {
        hpx::apply(boost::bind(
            &octree_client::tie_sibling_for_physical_boundary, this, _1, _2),
                target_f, target_sib);
        return;
    }

    hpx::apply<octree_server::tie_sibling_action>
        (gid_, target_f, target_sib);
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::set_child_sibling(
    child_index kid
  , face f
  , octree_client const& sib
    )
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::async<octree_server::set_child_sibling_action>
        (gid_, kid, f, sib).get();
}

void octree_client::set_child_sibling_push(
    child_index kid
  , face f
  , octree_client const& sib
    )
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::apply<octree_server::set_child_sibling_action>(gid_, kid, f, sib);
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::tie_child_sibling(
    child_index target_kid
  , face target_f
  , octree_client const& target_sib
    )
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::async<octree_server::set_child_sibling_action>
        (gid_, target_kid, target_f, target_sib).get();
}

void octree_client::tie_child_sibling_push(
    child_index target_kid
  , face target_f
  , octree_client const& target_sib
    )
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::apply<octree_server::set_child_sibling_action>
        (gid_, target_kid, target_f, target_sib);
}
    
///////////////////////////////////////////////////////////////////////////////
hpx::future<boost::array<octree_client, 6> > octree_client::get_siblings_async()
{
    ensure_real();
    return hpx::async<octree_server::get_siblings_action>(gid_);
}

///////////////////////////////////////////////////////////////////////////////
hpx::future<void> octree_client::receive_ghost_zones_async()
{
    ensure_real();
    return hpx::async<octree_server::receive_ghost_zones_action>(gid_);
}

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENT
vector3d<std::vector<double> > octree_client::interpolate(
    face f
    )
{
    return vector3d<std::vector<double> >();
}

// IMPLEMENT
vector3d<std::vector<double> > octree_client::mirror_or_outflow(
    face f
    )
{
    return vector3d<std::vector<double> >();
}

vector3d<std::vector<double> > octree_client::send_ghost_zone(
    face f
    )
{
    switch (kind_)
    {
        case real_boundary:
            return send_ghost_zone_async(f).get(); 
        case amr_boundary:
            return interpolate(f);
        case physical_boundary:
            return mirror_or_outflow(f); 
        default:
            break;
    }

    OCTOPUS_ASSERT(false);
    return vector3d<std::vector<double> >();
}

hpx::future<vector3d<std::vector<double> > > 
octree_client::send_ghost_zone_async(
    face f
    )
{
    switch (kind_)
    {
        case real_boundary:
            return hpx::async<octree_server::send_ghost_zone_action>(gid_, f);
        case amr_boundary:
            return hpx::async(
                boost::bind(&octree_client::interpolate, this, _1), f); 
        case physical_boundary:
            return hpx::async(
                boost::bind(&octree_client::mirror_or_outflow, this, _1), f);
        default:
            break;
    }

    OCTOPUS_ASSERT(false);
    return hpx::future<vector3d<std::vector<double> > >();
}

///////////////////////////////////////////////////////////////////////////////
hpx::future<void> octree_client::apply_async(
    hpx::util::function<void(octree_server&)> const& f
  , boost::uint64_t minimum_level
    )
{
    ensure_real();
    return hpx::async<octree_server::apply_action>(gid_, f, minimum_level);
}

void octree_client::apply_push(
    hpx::util::function<void(octree_server&)> const& f
  , boost::uint64_t minimum_level
    )
{
    ensure_real();
    hpx::apply<octree_server::apply_action>(gid_, f, minimum_level);
}

}

