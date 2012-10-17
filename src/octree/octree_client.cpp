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
#include <octopus/engine/engine_interface.hpp>
#include <octopus/operators/std_vector_arithmetic.hpp>
#include <octopus/math.hpp>

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
    ) const
{
    ensure_real();
    return hpx::async<octree_server::create_child_action>(gid_, kid);
}

void octree_client::set_sibling_for_amr_boundary(
    face f
  , octree_client const& sib 
  , octree_client const& parent
    ) const
{
    // IMPLEMENT
}

void octree_client::set_sibling_for_physical_boundary(
    face f
  , octree_client const& sib 
  , octree_client const& parent
    ) const
{
    // IMPLEMENT
}

///////////////////////////////////////////////////////////////////////////////
void octree_client::set_sibling(
    face f
  , octree_client const& sib
  , octree_client const& sib_parent
    ) const
{
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(f));

    if (amr_boundary == kind_)
    {
        set_sibling_for_amr_boundary(f, sib, sib_parent); 
        return;
    }

    else if (physical_boundary == kind_)
    {
        set_sibling_for_physical_boundary(f, sib, sib_parent); 
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
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > f,
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
        hpx::apply(boost::bind(
            &octree_client::set_sibling_for_physical_boundary, this,
                _1, _2, _3), f, sib, sib_parent);
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
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
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
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
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
    ) const
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
    ) const
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
    ) const
{
    ensure_real();
    OCTOPUS_ASSERT_FMT_MSG(out_of_bounds > target_f,
                           "invalid face, face(%1%)",
                           boost::uint16_t(target_f));
    hpx::apply<octree_server::set_child_sibling_action>
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
    vector3d<std::vector<double> > input =
        hpx::async<octree_server::send_ghost_zone_action>(gid_, f).get();

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().spatial_size;

    vector3d<std::vector<double> > output; 

    switch (f)
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

                        OCTOPUS_ASSERT(science().get_rho_from_state(m) > 0.0);
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

                        OCTOPUS_ASSERT(science().get_rho_from_state(m) > 0.0);
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

                        OCTOPUS_ASSERT(science().get_rho_from_state(m) > 0.0);
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

                        OCTOPUS_ASSERT(science().get_rho_from_state(m) > 0.0);
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

                        OCTOPUS_ASSERT(science().get_rho_from_state(m) > 0.0);
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

                        OCTOPUS_ASSERT(science().get_rho_from_state(m) > 0.0);
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

// IMPLEMENT
vector3d<std::vector<double> > octree_client::mirror_or_outflow(
    face f
    ) const
{
    return vector3d<std::vector<double> >();
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
    ) const
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
    ) const
{
    ensure_real();
    return hpx::async<octree_server::apply_action>(gid_, f, minimum_level);
}

void octree_client::apply_push(
    hpx::util::function<void(octree_server&)> const& f
  , boost::uint64_t minimum_level
    ) const
{
    ensure_real();
    hpx::apply<octree_server::apply_action>(gid_, f, minimum_level);
}

}

