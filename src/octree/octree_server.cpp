////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/config.hpp>
#include <hpx/async.hpp>
#include <hpx/lcos/future.hpp>
#include <hpx/lcos/future_wait.hpp>

#include <octopus/math.hpp>
#include <octopus/indexer2d.hpp>
#include <octopus/octree/octree_server.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/operators/boost_array_arithmetic.hpp>
#include <octopus/operators/std_vector_arithmetic.hpp>

// TODO: Add invariant checker for functions which should only be called during
// initialization. (I think I did this?).

namespace octopus
{

// TODO: Verify the size of parent_U and it's elements.
void octree_server::inject_state_from_parent(
    vector3d<std::vector<double> > const& parent_U 
    )
{ // {{{
    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().spatial_size;
    
    const indexer2d<2> indexer(bw, gnx - bw - 1, bw, gnx - bw - 1);

    mutex_type::scoped_lock l(mtx_);
  
    child_index c = get_child_index_locked(l);

    U_.resize(gnx);

    std::vector<double> s1(ss), s2(ss), s3(ss);

    for (boost::uint64_t index = 0; index <= indexer.maximum; ++index)
    {
        boost::uint64_t k = indexer.y(index);
        boost::uint64_t j = indexer.x(index);
        boost::uint64_t k0 = (bw + k) / 2 + c.z() * (gnx / 2 - bw);
        boost::uint64_t j0 = (bw + j) / 2 + c.y() * (gnx / 2 - bw);

        for ( boost::uint64_t i = bw, i0 = bw + c.x() * (gnx / 2 - bw)
            ; i < (gnx - bw)
            ; i += 2, ++i0)
        {
            std::vector<double> const& u = parent_U(i0, j0, k0);

            using namespace octopus::operators;

            s1 = minmod(parent_U(i0 + 1, j0, k0) - u
                      , u - parent_U(i0 - 1, j0, k0));

            s2 = minmod(parent_U(i0, j0 + 1, k0) - u
                      , u - parent_U(i0, j0 - 1, k0));

            s3 = minmod(parent_U(i0, j0, k0 + 1) - u
                      , u - parent_U(i0, j0, k0 - 1));

            // FIXME: The little DSEL makes for clean syntax, but I need to
            // check with Joel/Heller about how copy intensive this is.
            U_(i + 0, j + 0, k + 0) = u - (s1 + s2 + s3) * 0.25;
            U_(i + 1, j + 0, k + 0) = u + (s1 - s2 - s3) * 0.25;
            U_(i + 0, j + 1, k + 0) = u - (s1 - s2 + s3) * 0.25;
            U_(i + 1, j + 1, k + 0) = u + (s1 + s2 - s3) * 0.25;
            U_(i + 0, j + 0, k + 1) = u - (s1 + s2 - s3) * 0.25;
            U_(i + 1, j + 0, k + 1) = u + (s1 - s2 + s3) * 0.25;
            U_(i + 0, j + 1, k + 1) = u - (s1 - s2 - s3) * 0.25;
            U_(i + 1, j + 1, k + 1) = u + (s1 + s2 + s3) * 0.25;
        }
    }

    state_received_locked(l);
} // }}}

void octree_server::create_child(
    child_index kid
    )
{ // {{{
    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT_FMT_MSG(
        children_[kid] == hpx::naming::invalid_id,
        "child already exists, child(%1%)", kid);

    child_index x_sib = kid;
    child_index y_sib = kid;
    child_index z_sib = kid;

    OCTOPUS_TEST_IN_PLACE(x_sib == kid);
    OCTOPUS_TEST_IN_PLACE(y_sib == kid);
    OCTOPUS_TEST_IN_PLACE(z_sib == kid);

    // Exterior/interior is relative to the new child.
    face exterior_x_face = out_of_bounds; // f1 from original code.
    face interior_x_face = out_of_bounds; // f2 from original code.

    face exterior_y_face = out_of_bounds; // f1 from original code.
    face interior_y_face = out_of_bounds; // f2 from original code.

    face exterior_z_face = out_of_bounds; // f1 from original code.
    face interior_z_face = out_of_bounds; // f2 from original code.

    octree_init_data kid_init;

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().spatial_size;

    using namespace octopus::operators;

    kid_init.parent    = safe_reference(); 
    kid_init.level     = level_ + 1; 
    kid_init.location  = location_ * 2 + kid.array(); 
    kid_init.dx        = dx_ * 0.5;
    kid_init.time      = time_;
    kid_init.offset    = offset_ * 2 + bw + (kid.array() * (gnx - 2 * bw));
    kid_init.origin    = origin_;

    // Start creating the child. 
    hpx::future<hpx::id_type, hpx::naming::gid_type> kid_gid
        = create_octree_async(kid_init, U_);

    ///////////////////////////////////////////////////////////////////////////
    // X-axis. 
    if (0 == kid.x())
    {
        // The box that is in the (-1, 0, 0) direction (relative to this child)
        // is external, e.g. one of our siblings (or possibly an AMR/physics
        // boundary. The box that is in the (+1, 0, 0) direction is another
        // one of our siblings. 

        x_sib.set_x(1);

        OCTOPUS_TEST_IN_PLACE(x_sib.x() == 1);

        exterior_x_face = XL; // (-1, 0, 0)
        interior_x_face = XU; // (+1, 0, 0)
    }
    else
    {
        // The box that is in the (+1, 0, 0) direction (relative to this child)
        // is external, e.g. one of our siblings (or possibly an AMR/physics
        // boundary. The box that is in the (-1, 0, 0) direction is another
        // one of our siblings. 

        x_sib.set_x(0);

        OCTOPUS_TEST_IN_PLACE(x_sib.x() == 0);

        exterior_x_face = XU; // (+1, 0, 0)
        interior_x_face = XL; // (-1, 0, 0)
    }

    ///////////////////////////////////////////////////////////////////////////
    // Y-axis. 
    if (0 == kid.y())
    {
        // The box that is in the (0, -1, 0) direction (relative to this child)
        // is external, e.g. one of our siblings (or possibly an AMR/physics
        // boundary. The box that is in the (0, +1, 0) direction is another
        // one of our siblings. 

        y_sib.set_y(1);

        OCTOPUS_TEST_IN_PLACE(y_sib.y() == 1);

        exterior_y_face = YL; // (0, -1, 0)
        interior_y_face = YU; // (0, +1, 0)
    }
    else
    {
        // The box that is in the (0, +1, 0) direction (relative to this child)
        // is external, e.g. one of our siblings (or possibly an AMR/physics
        // boundary. The box that is in the (0, -1, 0) direction is another
        // one of our siblings. 

        y_sib.set_y(0);

        OCTOPUS_TEST_IN_PLACE(y_sib.y() == 0);

        exterior_y_face = YU; // (0, +1, 0)
        interior_y_face = YL; // (0, -1, 0)
    }

    ///////////////////////////////////////////////////////////////////////////
    // Z-axis. 
    if (0 == kid.z())
    {
        // The box that is in the (0, 0, -1) direction (relative to this child)
        // is external, e.g. one of our siblings (or possibly an AMR/physics
        // boundary. The box that is in the (0, 0, +1) direction is another
        // one of our siblings. 

        z_sib.set_z(1);

        OCTOPUS_TEST_IN_PLACE(z_sib.z() == 1);

        exterior_z_face = ZL; // (0, 0, -1)
        interior_z_face = ZU; // (0, 0, +1)
    }
    else
    {
        // The box that is in the (0, 0, +1) direction (relative to this child)
        // is external, e.g. one of our siblings (or possibly an AMR/physics
        // boundary. The box that is in the (0, 0, -1) direction is another
        // one of our siblings. 

        z_sib.set_z(0);

        OCTOPUS_TEST_IN_PLACE(z_sib.z() == 0);

        exterior_z_face = ZU; // (0, 0, +1)
        interior_z_face = ZL; // (0, 0, -1)
    }

    OCTOPUS_TEST_IN_PLACE(exterior_x_face != out_of_bounds);
    OCTOPUS_TEST_IN_PLACE(interior_x_face != out_of_bounds);
    OCTOPUS_TEST_IN_PLACE(exterior_y_face != out_of_bounds);
    OCTOPUS_TEST_IN_PLACE(interior_y_face != out_of_bounds);
    OCTOPUS_TEST_IN_PLACE(exterior_z_face != out_of_bounds);
    OCTOPUS_TEST_IN_PLACE(interior_z_face != out_of_bounds);

    // Now, we must wait for the child to be created.
    octree_client kid_client(kid_gid.get());

    OCTOPUS_ASSERT(kid_client != hpx::naming::invalid_id);

    ///////////////////////////////////////////////////////////////////////////
    // Create the interior "family" links.

    // Check if the interior X sibling of the new child exists.
    if (children_[x_sib] != hpx::naming::invalid_id)
        children_[x_sib].tie_sibling_push(exterior_x_face, kid_client);

    // Check if the interior Y sibling of the new child exists.
    if (children_[y_sib] != hpx::naming::invalid_id)
        children_[y_sib].tie_sibling_push(exterior_y_face, kid_client);

    // Check if the interior Z sibling of the new child exists.
    if (children_[z_sib] != hpx::naming::invalid_id)
        children_[z_sib].tie_sibling_push(exterior_z_face, kid_client);

    ///////////////////////////////////////////////////////////////////////////
    // Create the exterior "family" links.

    // Check if the exterior X uncle (get it? :D) of the new child exists.
    if (siblings_[exterior_x_face] != hpx::naming::invalid_id)
        siblings_[exterior_x_face].tie_child_sibling_push
            (x_sib, interior_x_face, kid_client); 

    // Check if the exterior Y uncle (get it? :D) of the new child exists.
    if (siblings_[exterior_y_face] != hpx::naming::invalid_id)
        siblings_[exterior_y_face].tie_child_sibling_push
            (y_sib, interior_y_face, kid_client); 

    // Check if the exterior Z uncle (get it? :D) of the new child exists.
    if (siblings_[exterior_z_face] != hpx::naming::invalid_id)
        siblings_[exterior_z_face].tie_child_sibling_push
            (z_sib, interior_z_face, kid_client); 
} // }}}

void octree_server::set_sibling(
    face f
  , octree_client const& sib
    )
{ // {{{
    OCTOPUS_ASSERT_FMT_MSG(
        out_of_bounds > f,
        "invalid face, face(%1%), sibling(%2%)",
        boost::uint16_t(f) % sib);

    {
        mutex_type::scoped_lock l(mtx_);
    
        OCTOPUS_ASSERT_FMT_MSG(
            siblings_[f] == hpx::naming::invalid_id,
            "sibling already exist, face(%1%), sibling(%2%)",
            boost::uint16_t(f) % sib);

        siblings_[f] = sib;  
        sibling_set_locked(l);
    }
} // }}}

void octree_server::tie_sibling(
    face target_f
  , octree_client target_sib
    )
{ // {{{
    // Locks.
    child_index target_kid = get_child_index();

    OCTOPUS_ASSERT_FMT_MSG(
        out_of_bounds > target_f,
        "invalid target face, target_kid(%1%), target_face(%2%), "
        "target_sibling(%3%)",
        target_kid % boost::uint16_t(target_f) % target_sib);

    child_index source_kid = target_kid;

    // Invert 
    switch (target_f)
    {
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        case XL: // source_kid = target_kid + (+1, 0, 0) 
        {
            OCTOPUS_ASSERT(target_kid.x() == 0);
            source_kid.set_x(1);
            break;
        } 
        case XU: // source_kid = target_kid + (-1, 0, 0) 
        {
            OCTOPUS_ASSERT(target_kid.x() == 1);
            source_kid.set_x(0);
            break;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        case YL: // source_kid = target_kid + (0, +1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.y() == 0);
            source_kid.set_y(1);
            break;
        } 
        case YU: // source_kid = target_kid + (0, -1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.y() == 1);
            source_kid.set_y(0);
            break;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        case ZL: // source_kid = target_kid + (0, +1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.z() == 0);
            source_kid.set_z(1);
            break;
        } 
        case ZU: // source_kid = target_kid + (0, -1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.z() == 1);
            source_kid.set_z(0);
            break;
        }

        default:
        {
            OCTOPUS_ASSERT_MSG(false, "source face shouldn't be out-of-bounds");
        }
    } 

    face source_f = invert(target_f);

    // Locks.
    set_sibling(target_f, target_sib);
    
    octree_client source_sib(get_gid());

    target_sib.set_sibling_push(source_f, source_sib);  
} // }}}

void octree_server::set_child_sibling(
    child_index kid
  , face f
  , octree_client const& sib
    )
{ // {{{
    OCTOPUS_ASSERT_FMT_MSG(
        out_of_bounds > f,
        "invalid face, kid(%1%), face(%2%), sibling(%3%)",
        kid % boost::uint16_t(f) % sib);

    // Make sure that we are initialized.
    initialized_.wait();

    {
        mutex_type::scoped_lock l(mtx_);

        OCTOPUS_ASSERT_FMT_MSG(
            children_[kid] != hpx::naming::invalid_id,
            "child does not exists, kid(%1%), face(%2%), sibling(%3%)",
            kid % boost::uint16_t(f) % sib);

        children_[kid].set_sibling_push(f, sib);
    }
} // }}}

void octree_server::tie_child_sibling(
    child_index target_kid
  , face target_f
  , octree_client target_sib
    )
{ // {{{
    OCTOPUS_ASSERT_FMT_MSG(
        out_of_bounds > target_f,
        "invalid target face, target_kid(%1%), target_face(%2%), "
        "target_sibling(%3%)",
        target_kid % boost::uint16_t(target_f) % target_sib);

    // Make sure that we are initialized.
    initialized_.wait();

    child_index source_kid = target_kid;

    // Invert 
    switch (target_f)
    {
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        case XL: // source_kid = target_kid + (+1, 0, 0) 
        {
            OCTOPUS_ASSERT(target_kid.x() == 0);
            source_kid.set_x(1);
            break;
        } 
        case XU: // source_kid = target_kid + (-1, 0, 0) 
        {
            OCTOPUS_ASSERT(target_kid.x() == 1);
            source_kid.set_x(0);
            break;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        case YL: // source_kid = target_kid + (0, +1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.y() == 0);
            source_kid.set_y(1);
            break;
        } 
        case YU: // source_kid = target_kid + (0, -1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.y() == 1);
            source_kid.set_y(0);
            break;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        case ZL: // source_kid = target_kid + (0, +1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.z() == 0);
            source_kid.set_z(1);
            break;
        } 
        case ZU: // source_kid = target_kid + (0, -1, 0) 
        {
            OCTOPUS_ASSERT(target_kid.z() == 1);
            source_kid.set_z(0);
            break;
        }

        default:
        {
            OCTOPUS_ASSERT_MSG(false, "source face shouldn't be out-of-bounds");
        }
    }; 

    face source_f = invert(target_f);

    octree_client source_sib;

    {
        mutex_type::scoped_lock l(mtx_);

        OCTOPUS_ASSERT_FMT_MSG(
            siblings_[source_f] != hpx::naming::invalid_id,
            "source face is a boundary, target_kid(%1%), target_face(%2%), "
            "source_kid(%3%), source_face(%4%), target_sibling(%5%)",
            target_kid % boost::uint16_t(target_f) %
            source_kid % boost::uint16_t(source_f) %
            target_sib);

        // Check if we're at a boundary.
        if (children_[target_kid] != hpx::naming::invalid_id)
        {
            children_[target_kid].set_sibling_push(target_f, target_sib);

            source_sib = children_[target_kid];
        }
    }

    // We established by assertion earlier that siblings_[source_f] is not NULL.
    // FIXME: If source_sib is NULL and that fact is implicitly known by the
    // caller, then this is non-optimal.
    siblings_[source_f].set_child_sibling_push
        (source_kid, source_f, children_[target_kid]);
} // }}}

boost::array<octree_client, 6> octree_server::get_siblings()
{ // {{{
    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);
    return siblings_;
} // }}}

void octree_server::inject_state_from_children()
{ // {{{ IMPLEMENT

} // }}}

///////////////////////////////////////////////////////////////////////////////
// Send/receive ghost zones
/// Pseudo code based on the original code (GS = GNX - 2 BW = dimensions of
/// the grid without ghostzones):
///
/// for i in [0, BW)
///     for j in [BW, GNX - BW)
///         for k in [BW, GNX - BW)
///             U(i, j, k) = sibling[XL].U(GNX - 2 * BW + i, j, k) 
///              
/// for i in [GNX - BW, GNX)
///     for j in [BW, GNX - BW)
///         for k in [BW, GNX - BW)
///             U(i, j, k) = sibling[XU].U(-GNX - 2 * BW + i, j, k) 
///
/// for i in [BW, GNX - BW)
///     for j in [0, BW)
///         for k in [BW, GNX - BW)
///             U(i, j, k) = sibling[YL].U(i, GNX - 2 * BW + j, k) 
///
/// for i in [BW, GNX - BW)
///     for j in [GNX - BW, GNX)
///         for k in [BW, GNX - BW)
///             U(i, j, k) = sibling[YU].U(i, -GNX - 2 * BW + j, k) 
///
/// for i in [BW, GNX - BW)
///     for j in [BW, GNX - BW)
///         for k in [0, BW)
///             U(i, j, k) = sibling[ZL].U(i, j, GNX - 2 * BW + k) 
///
/// for i in [BW, GNX - BW)
///     for j in [BW, GNX - BW)
///         for k in [GNX - BW, GNX)
///             U(i, j, k) = sibling[ZU].U(i, j, -GNX - 2 * BW + k) 

// Who ya gonna call? Ghostbusters!
vector3d<std::vector<double> > octree_server::send_ghost_zone(
    face f
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().spatial_size;

    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    switch (f)
    {
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        /// for i in [0, BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[XL].U(GNX - 2 * BW + i, j, k) 
        ///              
        case XL:
        {
            vector3d<std::vector<double> > zone
                (
                /* [0, BW) */         bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = 0; i < bw; ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i;
                        boost::uint64_t const jj = j - bw;
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(gnx - 2 * bw + i, j, k);
                    }

            return zone;
        } 

        /// for i in [GNX - BW, GNX)
        ///     for j in [BW, GNX - BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[XU].U(-GNX - 2 * BW + i, j, k) 
        case XU:
        {
            vector3d<std::vector<double> > zone
                (
                /* [GNX - BW, GNX) */ bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = gnx - bw; i < gnx; ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - (gnx - bw);
                        boost::uint64_t const jj = j - bw;
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(-gnx - 2 * bw + i, j, k);
                    }

            return zone;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        /// for i in [BW, GNX - BW)
        ///     for j in [0, BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[YL].U(i, GNX - 2 * BW + j, k) 
        ///
        case YL:
        {
            vector3d<std::vector<double> > zone
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [0, BW) */         bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = 0; j < bw; ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j;
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(i, gnx - 2 * bw + j, k);
                    }

            return zone;
        } 

        /// for i in [BW, GNX - BW)
        ///     for j in [GNX - BW, GNX)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[YU].U(i, -GNX - 2 * BW + j, k) 
        case YU:
        {
            vector3d<std::vector<double> > zone
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [GNX - BW, GNX) */ bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = gnx - bw; j < gnx; ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - (gnx - bw);
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(i, -gnx - 2 * bw + j, k);
                    }

            return zone;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        /// for i in [BW, GNX - BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [0, BW)
        ///             U(i, j, k) = sibling[ZL].U(i, j, GNX - 2 * BW + k) 
        case ZL:
        {
            vector3d<std::vector<double> > zone
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [0, BW) */         bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = 0; k < bw; ++k)
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k;

                        zone(ii, jj, kk) = U_(i, j, gnx - 2 * bw + k);
                    }

            return zone;
        } 

        /// for i in [BW, GNX - BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [GNX - BW, GNX)
        ///             U(i, j, k) = sibling[ZU].U(i, j, -GNX - 2 * BW + k) 
        case ZU:
        {
            vector3d<std::vector<double> > zone
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [GNX - BW, GNX) */ bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = gnx - bw; k < gnx; ++k)
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k - (gnx - bw);

                        zone(ii, jj, kk) = U_(i, j, -gnx - 2 * bw + k);
                    }

            return zone;
        }

        default:
        {
            OCTOPUS_ASSERT_MSG(false, "face shouldn't be out-of-bounds");
        }
    }; 

    // Unreachable.
    OCTOPUS_ASSERT(false);
    return vector3d<std::vector<double> >(); 
} // }}} 

void octree_server::integrate_ghost_zone(
    std::size_t f
  , vector3d<std::vector<double> > const& zone
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().spatial_size;

    // First, we need to re-acquire a lock on the mutex.
    mutex_type::scoped_lock l(mtx_);

    // The index of the futures in the vector is the face.
    switch (f)
    {
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        /// for i in [0, BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[XL].U(GNX - 2 * BW + i, j, k) 
        ///              
        case XL:
        {
            OCTOPUS_ASSERT(zone.x_length() == bw);           // [0, BW)
            OCTOPUS_ASSERT(zone.y_length() == gnx - 2 * bw); // [BW, GNX - BW)  
            OCTOPUS_ASSERT(zone.z_length() == gnx - 2 * bw); // [BW, GNX - BW)

            for (boost::uint64_t i = 0; i < bw; ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i;
                        boost::uint64_t const jj = j - bw;
                        boost::uint64_t const kk = k - bw; 

                        U_(i, j, k) = zone(ii, jj, kk);
                    }

            return;
        } 

        /// for i in [GNX - BW, GNX)
        ///     for j in [BW, GNX - BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[XU].U(-GNX - 2 * BW + i, j, k) 
        case XU:
        {
            OCTOPUS_ASSERT(zone.x_length() == bw);           // [GNX - BW, GNX)
            OCTOPUS_ASSERT(zone.y_length() == gnx - 2 * bw); // [BW, GNX - BW)  
            OCTOPUS_ASSERT(zone.z_length() == gnx - 2 * bw); // [BW, GNX - BW)

            for (boost::uint64_t i = gnx - bw; i < gnx; ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - (gnx - bw);
                        boost::uint64_t const jj = j - bw;
                        boost::uint64_t const kk = k - bw; 

                        U_(i, j, k) = zone(ii, jj, kk);
                    }

            return;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        /// for i in [BW, GNX - BW)
        ///     for j in [0, BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[YL].U(i, GNX - 2 * BW + j, k) 
        ///
        case YL:
        {
            OCTOPUS_ASSERT(zone.x_length() == gnx - 2 * bw); // [BW, GNX - BW)  
            OCTOPUS_ASSERT(zone.y_length() == bw);           // [0, BW)
            OCTOPUS_ASSERT(zone.z_length() == gnx - 2 * bw); // [BW, GNX - BW)

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = 0; j < bw; ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j;
                        boost::uint64_t const kk = k - bw; 

                        U_(i, j, k) = zone(ii, jj, kk);
                    }

            return;
        } 

        /// for i in [BW, GNX - BW)
        ///     for j in [GNX - BW, GNX)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[YU].U(i, -GNX - 2 * BW + j, k) 
        case YU:
        {
            OCTOPUS_ASSERT(zone.x_length() == gnx - 2 * bw); // [BW, GNX - BW)  
            OCTOPUS_ASSERT(zone.y_length() == bw);           // [GNX - BW, GNX)
            OCTOPUS_ASSERT(zone.z_length() == gnx - 2 * bw); // [BW, GNX - BW)

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = gnx - bw; j < gnx; ++j)
                    for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - (gnx - bw);
                        boost::uint64_t const kk = k - bw; 

                        U_(i, j, k) = zone(ii, jj, kk);
                    }

            return;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        /// for i in [BW, GNX - BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [0, BW)
        ///             U(i, j, k) = sibling[ZL].U(i, j, GNX - 2 * BW + k) 
        case ZL:
        {
            OCTOPUS_ASSERT(zone.x_length() == gnx - 2 * bw); // [BW, GNX - BW)  
            OCTOPUS_ASSERT(zone.y_length() == bw);           // [0, BW)
            OCTOPUS_ASSERT(zone.z_length() == gnx - 2 * bw); // [BW, GNX - BW)

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = 0; k < bw; ++k)
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k;

                        U_(i, j, k) = zone(ii, jj, kk);
                    }

            return;
        } 

        /// for i in [BW, GNX - BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [GNX - BW, GNX)
        ///             U(i, j, k) = sibling[ZU].U(i, j, -GNX - 2 * BW + k) 
        case ZU:
        {
            OCTOPUS_ASSERT(zone.x_length() == gnx - 2 * bw); // [BW, GNX - BW)  
            OCTOPUS_ASSERT(zone.y_length() == bw);           // [GNX - BW, GNX)
            OCTOPUS_ASSERT(zone.z_length() == gnx - 2 * bw); // [BW, GNX - BW)

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = gnx - bw; k < gnx; ++k)
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k - (gnx - bw);

                        U_(i, j, k) = zone(ii, jj, kk);
                    }

            return;
        }

        default:
        {
            OCTOPUS_ASSERT_MSG(false, "face shouldn't be out-of-bounds");
        }
    }; 
} // }}} 

void octree_server::receive_ghost_zones()
{ // {{{
    // Make sure that we are initialized.
    initialized_.wait();

    // REVIEW: It'd be deadlocky to call this function on an octree_server that
    // is it's own sibling. Is this ever possible? Can't safely check for it
    // without locking the mutex, so no point in adding an assert (it'd deadlock
    // before the assert fired).
    mutex_type::scoped_lock l(mtx_);

    // REVIEW: I believe doing this in parallel should be safe, because we are
    // reading from interior points and writing to ghost zone regions. So,
    // there should be no overlapping read/writes. I may be incorrect though.
    std::vector<hpx::future<void> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::naming::invalid_id != children_[i])
            recursion_is_parallelism.push_back
                (children_[i].receive_ghost_zones_async()); 

    OCTOPUS_ASSERT(hpx::naming::invalid_id != siblings_[XL]);
    OCTOPUS_ASSERT(hpx::naming::invalid_id != siblings_[XU]);
    OCTOPUS_ASSERT(hpx::naming::invalid_id != siblings_[YL]);
    OCTOPUS_ASSERT(hpx::naming::invalid_id != siblings_[YU]);
    OCTOPUS_ASSERT(hpx::naming::invalid_id != siblings_[ZL]);
    OCTOPUS_ASSERT(hpx::naming::invalid_id != siblings_[ZU]);

    // FIXME: Would be nice if hpx::wait took boost::arrays.
    std::vector<hpx::future<vector3d<std::vector<double> > > > ghostzones;
    ghostzones.reserve(6);

    // NOTE: send_ghost_zone_async does special client-side magic for physical
    // boundaries and AMR boundaries.
    ghostzones.push_back(siblings_[XL].send_ghost_zone_async(XL));
    ghostzones.push_back(siblings_[XU].send_ghost_zone_async(XU));
    ghostzones.push_back(siblings_[YL].send_ghost_zone_async(YL));
    ghostzones.push_back(siblings_[YU].send_ghost_zone_async(YU));
    ghostzones.push_back(siblings_[ZL].send_ghost_zone_async(ZL));
    ghostzones.push_back(siblings_[ZU].send_ghost_zone_async(ZU));

    {
        // Unlock the lock ... 
        hpx::util::unlock_the_lock<mutex_type::scoped_lock> ul(l);

        // ... start polling for our ghost zones ...
        // FIXME: Hartmut wants this reimplemented with hpx::wait_any.
        hpx::wait(ghostzones,
            boost::bind(&octree_server::integrate_ghost_zone, this, _1, _2));

        // ... and block while our children to receive their ghost zones.
        hpx::wait(recursion_is_parallelism); 
    }
} // }}}

void octree_server::apply(
    hpx::util::function<void(octree_server&)> const& f
  , boost::uint64_t minimum_level
    )
{ // {{{
    mutex_type::scoped_lock l(mtx_);

    std::vector<hpx::future<void> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::naming::invalid_id != children_[i])
            recursion_is_parallelism.push_back(
                children_[i].apply_async(f, minimum_level)); 

    // Invoke the function on ourselves.
    if (level_ >= minimum_level)
        f(*this);

    {
        // Unlock the lock ... 
        hpx::util::unlock_the_lock<mutex_type::scoped_lock> ul(l);

        // ... and block while our children to receive their ghost zones.
        hpx::wait(recursion_is_parallelism); 
    }
} // }}}

void octree_server::save_state()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::add_differentials(double dt, double beta)
{ // {{{ IMPLEMENT

} // }}}

void octree_server::clear_differentials()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::step(double dt)
{ // {{{
    OCTOPUS_ASSERT_MSG(0 != level_,
        "step may only be called on the root octree_server");

    OCTOPUS_ASSERT_MSG(0 < dt, "invalid timestep size");

    // U -> U0, recursively.
    save_state();

    // NOTE: I have no good place to put this, so I'm putting it here: we do
    // TVD RK3 (google is your friend).
    switch (config().runge_kutta_order)
    {
        case 1:
        {
            sub_step(dt, 1.0);
            inject_state_from_children();
            break;
        }

        case 2:
        {
            sub_step(dt, 1.0);
            inject_state_from_children();
            sub_step(dt, 0.5);
            inject_state_from_children();
            break; 
        }

        case 3:
        {
            sub_step(dt, 1.0);
            inject_state_from_children();
            sub_step(dt, 0.25);
            inject_state_from_children();
            sub_step(dt, 2.0 / 3.0);
            inject_state_from_children();
            break; 
        }

        default:
        {
            OCTOPUS_ASSERT_FMT_MSG(false,
                "runge-kutta order (%1%) is unsupported or invalid",
                config().runge_kutta_order);
        }
    };

    receive_ghost_zones();
    refine();

    ++step_;
    time_ += dt;
} // }}}

void octree_server::sub_step(double dt, double beta)
{ // {{{
    OCTOPUS_ASSERT_MSG(0 != level_,
        "sub_step may only be called on the root octree_server");

    receive_ghost_zones();

    clear_differentials();

    // FIXME: I think these could be computed in parallel, if they each had
    // their own F.

    compute_x_flux();
    adjust_x_flux();
    sum_x_differentials();

    compute_y_flux();
    adjust_y_flux();
    sum_y_differentials();

    compute_z_flux();
    adjust_z_flux();
    sum_z_differentials();

    add_differentials(dt, beta);
} // }}}

void octree_server::refine()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::compute_x_flux()
{ // {{{ IMPLEMENT
 
} // }}}

void octree_server::compute_y_flux()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::compute_z_flux()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::adjust_x_flux()
{ // {{{ IMPLEMENT
 
} // }}}

void octree_server::adjust_y_flux()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::adjust_z_flux()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::sum_x_differentials()
{ // {{{ IMPLEMENT
 
} // }}}

void octree_server::sum_y_differentials()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::sum_z_differentials()
{ // {{{ IMPLEMENT

} // }}}

}

