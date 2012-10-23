////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/lcos/future.hpp>
#include <hpx/lcos/future_wait.hpp>
#include <hpx/lcos/local/packaged_continuation.hpp>
#include <hpx/async.hpp>
#include <hpx/apply.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/runtime/components/runtime_support.hpp>

#include <octopus/octree/octree_server.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/operators/std_vector_arithmetic.hpp>
#include <octopus/operators/boost_array_arithmetic.hpp>
#include <octopus/math.hpp>

namespace octopus
{

///////////////////////////////////////////////////////////////////////////////
void octree_client::create_root(
    hpx::id_type const& locality
  , octree_init_data const& init
    ) 
{
    kind_ = real_boundary;

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
    kind_ = real_boundary;

    OCTOPUS_ASSERT_FMT_MSG(!(locality.get_msb() & 0xFF),
                           "target is not a locality, gid(%1%)",
                           locality);

    hpx::components::runtime_support rts(locality);
    gid_ = rts.create_component_async<octopus::octree_server>(init).get();
}

///////////////////////////////////////////////////////////////////////////////
hpx::future<void> octree_client::create_child_async(
    child_index kid
    ) const
{
    ensure_real();
    return hpx::async<octree_server::create_child_action>(gid_, kid);
}

// FIXME: Non-optimal, find a better way to get the offsets.
// P.S. A way that doesn't involve passing a billion parameters, e.g. something
// like *_init_data.
void octree_client::set_sibling_for_amr_boundary(
    face f
  , octree_client const& sib 
  , octree_client const& sib_parent
    ) const
{ // {{{
    gid_ = sib.gid_;

    // FIXME: Non-optimal.
    hpx::future<boost::array<boost::int64_t, 3> > sib_offset
        = sib.get_offset_async(); 

    // FIXME: Non-optimal.
    hpx::future<boost::array<boost::int64_t, 3> > sib_parent_offset
        = sib_parent.get_offset_async(); 

    face_ = f;

    boost::array<boost::int64_t, 3> v = { { 0, 0, 0 } };

    switch (f)
    {
        case XU:
            v[0] = -1;
            break;
        case XL:
            v[0] = 1;
            break;
        case YU:
            v[1] = -1;
            break;
        case YL:
            v[1] = 1;
            break;
        case ZU:
            v[2] = -1;
            break;
        case ZL:
            v[2] = 1;
            break;
        default:
            OCTOPUS_ASSERT(false);
            break;
    }

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    using namespace octopus::operators;

    v *= (gnx - 2 * bw);

    offset_ = sib_offset.get();
    offset_ += v;
    offset_ -= sib_parent_offset.get() * 2;
} // }}}

void octree_client::set_sibling_for_physical_boundary(
    face f
  , octree_client const& sib 
    ) const
{ // {{{
    gid_ = sib.gid_;

    switch (f)
    {
        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        case XU:
            reflect_ = config().x_reflect;
            direction_ = x_axis;
            return;
        case XL:
            reflect_ = false; 
            direction_ = x_axis;
            return;

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        case YU:
            reflect_ = config().y_reflect;
            direction_ = y_axis;
            return;
        case YL:
            reflect_ = false; 
            direction_ = y_axis;
            return;

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        case ZU:
            reflect_ = config().z_reflect;
            direction_ = z_axis;
            return;
        case ZL:
            reflect_ = false;
            direction_ = z_axis;
            return;
        default:
            break;
    }

    OCTOPUS_ASSERT(false);
} // }}}

///////////////////////////////////////////////////////////////////////////////
void octree_client::set_sibling(
    face f
  , octree_client const& sib
  , octree_client const& sib_parent
    ) const
{
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));

    if (amr_boundary == kind_)
    {
        set_sibling_for_amr_boundary(f, sib, sib_parent); 
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
  , octree_client const& sib_parent
    ) const
{
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));

    if (amr_boundary == kind_)
    {
        hpx::apply(boost::bind(
            &octree_client::set_sibling_for_amr_boundary, this,
                _1, _2, _3), f, sib, sib_parent);
        return;
    }

    else if (physical_boundary == kind_)
    {
        // This is guranteed to be a purely local operation, and is also
        // trivial, so we just do it directly.
        set_sibling_for_physical_boundary(f, sib); 
        return;
    }

    hpx::apply<octree_server::set_sibling_action>(gid_, f, sib);
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::tie_sibling(
    face target_f
  , octree_client const& target_sib
  , octree_client const& target_sib_parent
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::async<octree_server::tie_sibling_action>
        (gid_, target_f, target_sib, target_sib_parent).get();
}

void octree_client::tie_sibling_push(
    face target_f
  , octree_client const& target_sib
  , octree_client const& target_sib_parent
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::apply<octree_server::tie_sibling_action>
        (gid_, target_f, target_sib, target_sib_parent);
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::set_child_sibling(
    child_index kid
  , face f
  , octree_client const& sib
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::async<octree_server::set_child_sibling_action>
        (gid_, kid, f, sib).get();
}

void octree_client::set_child_sibling_push(
    child_index kid
  , face f
  , octree_client const& sib
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));
    hpx::apply<octree_server::set_child_sibling_action>(gid_, kid, f, sib);
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::tie_child_sibling(
    child_index target_kid
  , face target_f
  , octree_client const& target_sib
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::async<octree_server::tie_child_sibling_action>
        (gid_, target_kid, target_f, target_sib).get();
}

void octree_client::tie_child_sibling_push(
    child_index target_kid
  , face target_f
  , octree_client const& target_sib
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(invalid_face > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::apply<octree_server::tie_child_sibling_action>
        (gid_, target_kid, target_f, target_sib);
}
    
///////////////////////////////////////////////////////////////////////////////
hpx::future<boost::array<octree_client, 6> >
octree_client::get_siblings_async() const
{
    ensure_real();
    return hpx::async<octree_server::get_siblings_action>(gid_);
}

///////////////////////////////////////////////////////////////////////////////
hpx::future<boost::array<boost::int64_t, 3> >
octree_client::get_offset_async() const
{
    ensure_real();
    return hpx::async<octree_server::get_offset_action>(gid_);
}

///////////////////////////////////////////////////////////////////////////////
hpx::future<void> octree_client::receive_ghost_zones_async() const
{
    ensure_real();
    return hpx::async<octree_server::receive_ghost_zones_action>(gid_);
}

///////////////////////////////////////////////////////////////////////////////
// FIXME: This could use .when() continuations to be more asynchronous, maybe.
// FIXME: Only get the data we need from the parent.
// FIXME: Interpolate in place? Is this possible?
vector3d<std::vector<double> > octree_client::interpolate(
    face f
    ) const
{ // {{{
    // set_sibling(f), f is the direction of the caller relative to the sibling
    // (the sibling == us). send_ghost_zone, f is our direction relative to the
    // the caller. REVIEW: I think.
    OCTOPUS_ASSERT(face_ == invert(f));

    vector3d<std::vector<double> > input =
        hpx::async<octree_server::send_ghost_zone_action>(gid_, face_).get();

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    vector3d<std::vector<double> > output; 

    switch (face_)
    {
        
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        case XL:
        {
            output.resize
                (
                /* [0, BW) */         bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = 0; i < bw; ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        using namespace octopus::operators;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (output). 
                        boost::uint64_t const i_out = i;
                        boost::uint64_t const j_out = j - bw;
                        boost::uint64_t const k_out = k - bw; 

                        ///////////////////////////////////////////////////////
                        bool const i0 = (offset_[0] + i) % 2;

                        boost::uint64_t const i1 = (offset_[0] + i) / 2;
                        boost::uint64_t const j1 = (offset_[1] + j) / 2;
                        boost::uint64_t const k1 = (offset_[2] + k) / 2;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (input). 
                        boost::uint64_t const i_in = i1;
                        boost::uint64_t const j_in = j1 - bw;
                        boost::uint64_t const k_in = k1 - bw; 

                        std::vector<double> const& u = input(i_in, j_in, k_in); 
                        std::vector<double>& m = output(i_out, j_out, k_out); 

                        m = minmod(input(i_in + 1, j_in, k_in) - u
                                 , u - input(i_in - 1, j_in, k_in));

                        if (1 == i0)
                            m = -m;

                        m -= m * 0.25;

                        // FIXME: This is too specific to Dominic/Zach's
                        // code, move this into the science table.
                        //OCTOPUS_ASSERT(science().rho(m) > 0.0);
                    }

            break;
        } 

        case XU:
        {
            output.resize
                (
                /* [GNX - BW, GNX) */ bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = gnx - bw; i < gnx; ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        using namespace octopus::operators;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (output). 
                        boost::uint64_t const i_out = i - (gnx - bw);
                        boost::uint64_t const j_out = j - bw;
                        boost::uint64_t const k_out = k - bw; 

                        ///////////////////////////////////////////////////////
                        bool const i0 = (offset_[0] + i) % 2;

                        boost::uint64_t const i1 = (offset_[0] + i) / 2;
                        boost::uint64_t const j1 = (offset_[1] + j) / 2;
                        boost::uint64_t const k1 = (offset_[2] + k) / 2;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (input). 
                        boost::uint64_t const i_in = i1;
                        boost::uint64_t const j_in = j1 - bw;
                        boost::uint64_t const k_in = k1 - bw; 

                        std::vector<double> const& u = input(i_in, j_in, k_in); 
                        std::vector<double>& m = output(i_out, j_out, k_out); 

                        m = minmod(input(i_in + 1, j_in, k_in) - u
                                 , u - input(i_in - 1, j_in, k_in));

                        if (1 == i0)
                            m = -m;

                        m -= m * 0.25;

                        // FIXME: This is too specific to Dominic/Zach's
                        // code, move this into the science table.
                        //OCTOPUS_ASSERT(science().rho(m) > 0.0);
                    }

            break;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        case YL:
        {
            output.resize
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [0, BW) */         bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = 0; j < bw; ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        using namespace octopus::operators;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (output). 
                        boost::uint64_t const i_out = i - bw;
                        boost::uint64_t const j_out = j;
                        boost::uint64_t const k_out = k - bw; 

                        ///////////////////////////////////////////////////////
                        bool const j0 = (offset_[1] + j) % 2;

                        boost::uint64_t const i1 = (offset_[0] + i) / 2;
                        boost::uint64_t const j1 = (offset_[1] + j) / 2;
                        boost::uint64_t const k1 = (offset_[2] + k) / 2;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (input). 
                        boost::uint64_t const i_in = i1 - bw;
                        boost::uint64_t const j_in = j1;
                        boost::uint64_t const k_in = k1 - bw; 

                        std::vector<double> const& u = input(i_in, j_in, k_in); 
                        std::vector<double>& m = output(i_out, j_out, k_out); 

                        m = minmod(input(i_in, j_in + 1, k_in) - u
                                 , u - input(i_in, j_in - 1, k_in));

                        if (1 == j0)
                            m = -m;

                        m -= m * 0.25;

                        // FIXME: This is too specific to Dominic/Zach's
                        // code, move this into the science table.
                        //OCTOPUS_ASSERT(science().rho(m) > 0.0);
                    }

            break;
        } 

        case YU:
        {
            output.resize
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [GNX - BW, GNX) */ bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = gnx - bw; j < gnx; ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        using namespace octopus::operators;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (output). 
                        boost::uint64_t const i_out = i - bw;
                        boost::uint64_t const j_out = j - (gnx - bw);
                        boost::uint64_t const k_out = k - bw; 

                        ///////////////////////////////////////////////////////
                        bool const j0 = (offset_[1] + j) % 2;

                        boost::uint64_t const i1 = (offset_[0] + i) / 2;
                        boost::uint64_t const j1 = (offset_[1] + j) / 2;
                        boost::uint64_t const k1 = (offset_[2] + k) / 2;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (input). 
                        boost::uint64_t const i_in = i1 - bw;
                        boost::uint64_t const j_in = j1;
                        boost::uint64_t const k_in = k1 - bw; 

                        std::vector<double> const& u = input(i_in, j_in, k_in); 
                        std::vector<double>& m = output(i_out, j_out, k_out); 

                        m = minmod(input(i_in, j_in + 1, k_in) - u
                                 , u - input(i_in, j_in - 1, k_in));

                        if (1 == j0)
                            m = -m;

                        m -= m * 0.25;

                        // FIXME: This is too specific to Dominic/Zach's
                        // code, move this into the science table.
                        //OCTOPUS_ASSERT(science().rho(m) > 0.0);
                    }

            break;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        case ZL:
        {
            output.resize
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [0, BW) */         bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = 0; k < bw; ++k)
                    {
                        using namespace octopus::operators;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (output). 
                        boost::uint64_t const i_out = i - bw;
                        boost::uint64_t const j_out = j - bw; 
                        boost::uint64_t const k_out = k;

                        ///////////////////////////////////////////////////////
                        bool const k0 = (offset_[1] + k) % 2;

                        boost::uint64_t const i1 = (offset_[0] + i) / 2;
                        boost::uint64_t const j1 = (offset_[1] + j) / 2;
                        boost::uint64_t const k1 = (offset_[2] + k) / 2;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (input). 
                        boost::uint64_t const i_in = i1 - bw;
                        boost::uint64_t const j_in = j1 - bw; 
                        boost::uint64_t const k_in = k1;

                        std::vector<double> const& u = input(i_in, j_in, k_in); 
                        std::vector<double>& m = output(i_out, j_out, k_out); 

                        m = minmod(input(i_in, j_in, k_in + 1) - u
                                 , u - input(i_in, j_in, k_in - 1));

                        if (1 == k0)
                            m = -m;

                        m -= m * 0.25;

                        // FIXME: This is too specific to Dominic/Zach's
                        // code, move this into the science table.
                        //OCTOPUS_ASSERT(science().rho(m) > 0.0);
                    }

            break;
        } 

        case ZU:
        {
            output.resize
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [GNX - BW, GNX) */ bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = gnx - bw; k < gnx; ++k)
                    {
                        using namespace octopus::operators;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (output). 
                        boost::uint64_t const i_out = i - bw;
                        boost::uint64_t const j_out = j - bw; 
                        boost::uint64_t const k_out = k - (gnx - bw);

                        ///////////////////////////////////////////////////////
                        bool const k0 = (offset_[2] + k) % 2;

                        boost::uint64_t const i1 = (offset_[0] + i) / 2;
                        boost::uint64_t const j1 = (offset_[1] + j) / 2;
                        boost::uint64_t const k1 = (offset_[2] + k) / 2;

                        ///////////////////////////////////////////////////////
                        // Adjusted indices (input). 
                        boost::uint64_t const i_in = i1 - bw;
                        boost::uint64_t const j_in = j1 - bw; 
                        boost::uint64_t const k_in = k1;

                        std::vector<double> const& u = input(i_in, j_in, k_in); 
                        std::vector<double>& m = output(i_out, j_out, k_out); 

                        m = minmod(input(i_in, j_in, k_in + 1) - u
                                 , u - input(i_in, j_in, k_in - 1));

                        if (1 == k0)
                            m = -m;

                        m -= m * 0.25;

                        // FIXME: This is too specific to Dominic/Zach's
                        // code, move this into the science table.
                        //OCTOPUS_ASSERT(science().rho(m) > 0.0);
                    }

            break;
        }

        default:
        {
            OCTOPUS_ASSERT_MSG(false, "face shouldn't be out-of-bounds");
        }
    }; 

    return output; 
} // }}}

hpx::future<vector3d<std::vector<double> > > octree_client::interpolate_async(
    face f
    ) const
{
    return hpx::async(boost::bind(&octree_client::interpolate, this, _1), f); 
}

vector3d<std::vector<double> > octree_client::map(
    face f
    ) const
{
    // set_sibling(f), f is the direction of the caller relative to the sibling
    // (the sibling == us). send_ghost_zone, f is our direction relative to the
    // the caller. REVIEW: I think.
    OCTOPUS_ASSERT(face_ == invert(f));
    return hpx::async<octree_server::send_mapped_ghost_zone_action>
        (gid_, face_, reflect_, direction_).get();
}

hpx::future<vector3d<std::vector<double> > > octree_client::map_async(
    face f
    ) const
{
    // set_sibling(f), f is the direction of the caller relative to the sibling
    // (the sibling == us). send_ghost_zone, f is our direction relative to the
    // the caller. REVIEW: I think.
    OCTOPUS_ASSERT(face_ == invert(f));
    return hpx::async<octree_server::send_mapped_ghost_zone_action>
        (gid_, face_, reflect_, direction_);
}

vector3d<std::vector<double> > octree_client::send_ghost_zone(
    face f
    ) const
{
    switch (kind_)
    {
        case real_boundary:
            return send_ghost_zone_async(f).get(); 
        case amr_boundary:
            return interpolate(f);
        case physical_boundary:
            return map(f); 
        default:
            break;
    }

    OCTOPUS_ASSERT(false);
    return vector3d<std::vector<double> >();
}

hpx::future<vector3d<std::vector<double> > > 
octree_client::send_ghost_zone_async(
    face f
    ) const
{
    switch (kind_)
    {
        case real_boundary:
            return hpx::async<octree_server::send_ghost_zone_action>(gid_, f);
        case amr_boundary:
            return interpolate_async(f);
        case physical_boundary:
            return map_async(f); 
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
    ) const
{
    ensure_real();
    return hpx::async<octree_server::apply_action>
        (gid_, f, minimum_level);
}

void octree_client::apply_push(
    hpx::util::function<void(octree_server&)> const& f
  , boost::uint64_t minimum_level
    ) const
{
    ensure_real();
    hpx::apply<octree_server::apply_action>(gid_, f, minimum_level);
}

hpx::future<void> octree_client::step_async(double dt) const
{
    ensure_real();
    return hpx::async<octree_server::step_action>(gid_, dt);
}

void octree_client::step_push(double dt) const
{
    ensure_real();
    hpx::apply<octree_server::step_action>(gid_, dt);
}

void octree_client::step_to_time_push(double dt, double until) const
{
    ensure_real();
    hpx::apply<octree_server::step_to_time_action>(gid_, dt, until);
}

void begin_io_epoch(boost::uint64_t step)
{
    science().output.open(step);
}

void end_io_epoch()
{
    science().output.close();
}

}

HPX_PLAIN_ACTION(octopus::begin_io_epoch, begin_io_epoch_action);
HPX_PLAIN_ACTION(octopus::end_io_epoch, end_io_epoch_action);

namespace octopus
{

struct end_io_epoch_continuation
{
    typedef void result_type;

    hpx::future<void> f_;

    end_io_epoch_continuation(hpx::future<void> const& f)
      : f_(f)
    {}

    result_type operator()(hpx::future<void> res) const
    {
        std::vector<hpx::id_type> const& targets = localities();

        std::vector<hpx::future<void> > futures;
        futures.reserve(targets.size());

        end_io_epoch_action act;

        for (boost::uint64_t i = 0; i < targets.size(); ++i)
            futures.push_back(hpx::async(act, targets[i]));

        hpx::wait(futures);
    }
};

hpx::future<void> octree_client::output_async() const
{
    ensure_real();

    ///////////////////////////////////////////////////////////////////////////
    // Begin I/O epoch.
    // FIXME: Sub-optimal.
    hpx::future<boost::uint64_t> step_ 
        = hpx::async<octree_server::get_step_action>(gid_);

    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    begin_io_epoch_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
        futures.push_back(hpx::async(act, targets[i], step_.get()));

    hpx::wait(futures);

    ///////////////////////////////////////////////////////////////////////////
    // FIXME: This is broken.
    //return hpx::async<octree_server::output_action>(gid_).when
    //    (end_io_epoch_continuation());

    hpx::future<void> f = hpx::async<octree_server::output_action>(gid_);
    return f.when(end_io_epoch_continuation(f));
}

}

