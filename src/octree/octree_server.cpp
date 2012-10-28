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
#include <hpx/lcos/wait_all.hpp>

#include <octopus/math.hpp>
#include <octopus/indexer2d.hpp>
#include <octopus/octree/octree_server.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/operators/boost_array_arithmetic.hpp>
#include <octopus/operators/std_vector_arithmetic.hpp>

// TODO: Verify the size of parent_U and it's elements when initialization is
// complete.

// NOTE (wash): Is it necessary to solve coarser regions of the grid that are
// being solved at a finer level? I know this is necessary for the multigrid
// solver for the Poisson in the original binary code, but do we need it for
// the non-self gravating code?

namespace octopus
{

void octree_server::inject_state_from_parent(
    vector3d<std::vector<double> > const& parent_U 
    )
{ // {{{
    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;
    
    indexer2d<2> const indexer(bw, gnx - bw - 1, bw, gnx - bw - 1);

    mutex_type::scoped_lock l(mtx_);
  
    child_index c = get_child_index_locked(l);

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
            // check with Joel Falcou/Heller about how copy intensive this is.
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

void octree_server::prepare_queues()
{
    OCTOPUS_ASSERT(1 <= config().runge_kutta_order);
    OCTOPUS_ASSERT(3 >= config().runge_kutta_order);

    // NOTE: See the math in the header (right before the declaration of
    // ghost_zone_queue_) to see where these numbers come from. 
    ghost_zone_queue_.resize(
        config().runge_kutta_order + 1
      , sibling_dependencies());

    if (config().max_refinement_level == level_)
        return;

    children_state_queue_.resize(
        config().runge_kutta_order
      , children_dependencies());

    adjust_flux_queue_.resize(
        config().runge_kutta_order
      , children_dependencies()
        );
}


/// \brief Construct a root node. 
octree_server::octree_server(
    back_pointer_type back_ptr
  , octree_init_data const& init
    )
  : base_type(back_ptr)
  , initialized_()
  , mtx_()
  , siblings_set_(6)
  , state_received_(true)
  , future_self_(init.future_self)
  , past_self_(init.past_self)
  , ghost_zone_queue_()
  , children_state_queue_()
  , adjust_flux_queue_()
  , parent_(init.parent)
  , siblings_()
  , level_(init.level)
  , location_(init.location)
  , dx_(init.dx)
  , dx0_(science().initial_spacestep())
  , time_(init.time)
  , offset_(init.offset)
  , origin_(init.origin)
  , step_(0)
  , U_(config().grid_node_length, std::vector<double>(science().state_size))
  , U0_()
  , FX_(config().grid_node_length, std::vector<double>(science().state_size))
  , FY_(config().grid_node_length, std::vector<double>(science().state_size))
  , FZ_(config().grid_node_length, std::vector<double>(science().state_size))
  , FO_(science().state_size)
  , FO0_()
  , D_(config().grid_node_length, std::vector<double>(science().state_size))
  , DFO_(science().state_size)
{
    OCTOPUS_ASSERT(parent_ == hpx::invalid_id);

    prepare_queues();

    initialized_.set();

    for (face i = XL; i < invalid_face; i = face(boost::uint8_t(i + 1)))
    {
        siblings_[i] = octree_client(physical_boundary); 
        siblings_[i].set_sibling_for_physical_boundary(i, client_from_this());
    } 
}

// FIXME: Non-optimal, inject_state_from_parent should be called with
// fire-and-forget semantics. However, the lifetime of parent_U then
// becomes an issue. Because of this, U_ may need to be made into a
// shared_ptr.
/// \brief Construct a child node.
octree_server::octree_server(
    back_pointer_type back_ptr
  , octree_init_data const& init
  , vector3d<std::vector<double> > const& parent_U
    )
  : base_type(back_ptr)
  , initialized_()
  , mtx_()
  , siblings_set_(0)
  , state_received_(false)
  , future_self_(init.future_self)
  , past_self_(init.past_self)
  , ghost_zone_queue_()
  , children_state_queue_()
  , adjust_flux_queue_()
  , parent_(init.parent)
  , siblings_()
  , level_(init.level)
  , location_(init.location)
  , dx_(init.dx)
  , dx0_(science().initial_spacestep())
  , time_(init.time)
  , offset_(init.offset)
  , origin_(init.origin)
  , step_(init.step)
  , U_(config().grid_node_length, std::vector<double>(science().state_size))
  , U0_()
  , FX_(config().grid_node_length, std::vector<double>(science().state_size))
  , FY_(config().grid_node_length, std::vector<double>(science().state_size))
  , FZ_(config().grid_node_length, std::vector<double>(science().state_size))
  , FO_(science().state_size)
  , FO0_()
  , D_(config().grid_node_length, std::vector<double>(science().state_size))
  , DFO_(science().state_size)
{
    // Make sure our parent reference is not reference counted.
    OCTOPUS_ASSERT_MSG(
        init.parent.get_management_type() == hpx::id_type::unmanaged,
        "reference cycle detected in child");

    prepare_queues();

    inject_state_from_parent(parent_U);
}

// NOTE: Should be thread-safe, offset_ and origin_ are only read, and never
// written to.
double octree_server::x_face(boost::uint64_t i) const
{
    boost::uint64_t const bw = science().ghost_zone_width;
    double const grid_dim = config().spatial_domain;

    using namespace octopus::operators;

    return double(offset_[0] + i) * dx_ - grid_dim - bw * dx0_ - origin_[0];
}


// NOTE: Should be thread-safe, offset_ and origin_ are only read, and never
// written to.
double octree_server::y_face(boost::uint64_t i) const
{
    boost::uint64_t const bw = science().ghost_zone_width;
    double const grid_dim = config().spatial_domain;

    using namespace octopus::operators;

    return double(offset_[1] + i) * dx_ - grid_dim - bw * dx0_ - origin_[1];
}

// NOTE: Should be thread-safe, offset_ and origin_ are only read, and never
// written to.
double octree_server::z_face(boost::uint64_t i) const
{
    boost::uint64_t const bw = science().ghost_zone_width;
    double const grid_dim = config().spatial_domain;

    using namespace octopus::operators;

    if (config().reflect_on_z)
        return double(offset_[2] + i) * dx_ - bw * dx0_ - origin_[2];
    else
        return double(offset_[2] + i) * dx_ - grid_dim - bw * dx0_;
}

void octree_server::create_child(
    child_index kid
    )
{ // {{{
    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT_FMT_MSG(
        children_[kid] == hpx::invalid_id,
        "child already exists, child(%1%)", kid);

    child_index x_sib = kid;
    child_index y_sib = kid;
    child_index z_sib = kid;

    OCTOPUS_ASSERT(x_sib == kid);
    OCTOPUS_ASSERT(y_sib == kid);
    OCTOPUS_ASSERT(z_sib == kid);

    // Exterior/interior is relative to the new child.
    face exterior_x_face = invalid_face; // f1 from original code.
    face interior_x_face = invalid_face; // f2 from original code.

    face exterior_y_face = invalid_face; // f1 from original code.
    face interior_y_face = invalid_face; // f2 from original code.

    face exterior_z_face = invalid_face; // f1 from original code.
    face interior_z_face = invalid_face; // f2 from original code.

    octree_init_data kid_init;

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    using namespace octopus::operators;

    kid_init.parent         = reference_from_this(); 
    kid_init.level          = level_ + 1; 
    kid_init.location       = location_ * 2 + kid.array(); 
    kid_init.dx             = dx_ * 0.5;
    kid_init.time           = time_;
    kid_init.offset         = offset_ * 2 + bw + (kid.array() * (gnx - 2 * bw));
    kid_init.origin         = origin_;
    kid_init.step           = step_;

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

        OCTOPUS_ASSERT(x_sib.x() == 1);

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

        OCTOPUS_ASSERT(x_sib.x() == 0);

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

        OCTOPUS_ASSERT(y_sib.y() == 1);

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

        OCTOPUS_ASSERT(y_sib.y() == 0);

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

        OCTOPUS_ASSERT(z_sib.z() == 1);

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

        OCTOPUS_ASSERT(z_sib.z() == 0);

        exterior_z_face = ZU; // (0, 0, +1)
        interior_z_face = ZL; // (0, 0, -1)
    }

    OCTOPUS_ASSERT(exterior_x_face != invalid_face);
    OCTOPUS_ASSERT(interior_x_face != invalid_face);
    OCTOPUS_ASSERT(exterior_y_face != invalid_face);
    OCTOPUS_ASSERT(interior_y_face != invalid_face);
    OCTOPUS_ASSERT(exterior_z_face != invalid_face);
    OCTOPUS_ASSERT(interior_z_face != invalid_face);

    // Now, we must wait for the child to be created.
    octree_client kid_client(kid_gid.get());

    OCTOPUS_ASSERT(kid_client != hpx::invalid_id);

    ///////////////////////////////////////////////////////////////////////////
    // Create the interior "family" links.

    // These links will be set later if the child doesn't exist yet. We assume
    // that if one child is being created, then all 7 other children will be
    // created. Do not violate this assumption.

    octree_client this_ = client_from_this();

    // Check if the interior X sibling of the new child exists.
    if (children_[x_sib] != hpx::invalid_id)
        children_[x_sib].tie_sibling_push(exterior_x_face, kid_client, this_);

    // Check if the interior Y sibling of the new child exists.
    if (children_[y_sib] != hpx::invalid_id)
        children_[y_sib].tie_sibling_push(exterior_y_face, kid_client, this_);

    // Check if the interior Z sibling of the new child exists.
    if (children_[z_sib] != hpx::invalid_id)
        children_[z_sib].tie_sibling_push(exterior_z_face, kid_client, this_);

    ///////////////////////////////////////////////////////////////////////////
    // Create the exterior "family" links.

    // These links must exist. They may be non-real (e.g. boundaries), but they
    // must exist.

    // Check if the exterior X uncle (get it? :D) of the new child exists.
    OCTOPUS_ASSERT(siblings_[exterior_x_face] != hpx::invalid_id);

    switch (siblings_[exterior_x_face].kind())
    {
        case real_boundary:
        {
            siblings_[exterior_x_face].tie_child_sibling_push
                (x_sib, interior_x_face, kid_client);
            break;
        }

        case amr_boundary:
        {
            octree_client bound(amr_boundary); 

            kid_client.tie_sibling_push(interior_x_face, bound,
                client_from_this());

            break;
        }

        case physical_boundary:
        {
            octree_client bound(physical_boundary); 

            kid_client.tie_sibling_push(interior_x_face, bound,
                client_from_this());

            break;
        }

        default:
        {
            OCTOPUS_ASSERT(false);
            break;
        }
    };

    // Check if the exterior Y uncle (get it? :D) of the new child exists.
    OCTOPUS_ASSERT(siblings_[exterior_y_face] != hpx::invalid_id);

    switch (siblings_[exterior_y_face].kind())
    {
        case real_boundary:
        {
            siblings_[exterior_y_face].tie_child_sibling_push
                (y_sib, interior_y_face, kid_client);
            break;
        }

        case amr_boundary:
        {
            octree_client bound(amr_boundary); 

            kid_client.tie_sibling_push(interior_y_face, bound,
                client_from_this());

            break;
        }

        case physical_boundary:
        {
            octree_client bound(physical_boundary); 

            kid_client.tie_sibling_push(interior_y_face, bound,
                client_from_this());

            break;
        }

        default:
        {
            OCTOPUS_ASSERT(false);
            break;
        }
    };

    // Check if the exterior Z uncle (get it? :D) of the new child exists.
    OCTOPUS_ASSERT(siblings_[exterior_z_face] != hpx::invalid_id);

    switch (siblings_[exterior_z_face].kind())
    {
        case real_boundary:
        {
            siblings_[exterior_z_face].tie_child_sibling_push
                (z_sib, interior_z_face, kid_client);
            break;
        }

        case amr_boundary:
        {
            octree_client bound(amr_boundary); 

            kid_client.tie_sibling_push(interior_z_face, bound,
                client_from_this());

            break;
        }

        case physical_boundary:
        {
            octree_client bound(physical_boundary); 

            kid_client.tie_sibling_push(interior_z_face, bound,
                client_from_this());

            break;
        }

        default:
        {
            OCTOPUS_ASSERT(false);
            break;
        }
    };
} // }}}

void octree_server::set_sibling(
    face f
  , octree_client const& sib
    )
{ // {{{
    OCTOPUS_ASSERT_FMT_MSG(
        invalid_face > f,
        "invalid face, face(%1%), sibling(%2%)",
        boost::uint16_t(f) % sib);

    {
        mutex_type::scoped_lock l(mtx_);
    
        OCTOPUS_ASSERT_FMT_MSG(
            siblings_[f] == hpx::invalid_id,
            "sibling already exist, face(%1%), sibling(%2%)",
            boost::uint16_t(f) % sib);

        siblings_[f] = sib;  
        sibling_set_locked(l);
    }
} // }}}

void octree_server::tie_sibling(
    face target_f
  , octree_client const& target_sib
  , octree_client const& target_sib_parent
    )
{ // {{{
    // Locks.
    child_index target_kid = get_child_index();

    OCTOPUS_ASSERT_FMT_MSG(
        invalid_face > target_f,
        "invalid target face, target_kid(%1%), target_face(%2%), "
        "target_sibling(%3%)",
        target_kid % boost::uint16_t(target_f) % target_sib);

    child_index source_kid = invert(target_f, target_kid);

    face source_f = invert(target_f);

    // Locks.
    set_sibling(target_f, target_sib);
    
    octree_client source_sib(get_gid());

    target_sib.set_sibling_push(source_f, source_sib, target_sib_parent);  
} // }}}

void octree_server::set_child_sibling(
    child_index kid
  , face f
  , octree_client const& sib
    )
{ // {{{
    OCTOPUS_ASSERT_FMT_MSG(
        invalid_face > f,
        "invalid face, kid(%1%), face(%2%), sibling(%3%)",
        kid % boost::uint16_t(f) % sib);

    // Make sure that we are initialized.
    initialized_.wait();

    {
        mutex_type::scoped_lock l(mtx_);

        OCTOPUS_ASSERT_FMT_MSG(
            children_[kid] != hpx::invalid_id,
            "child does not exists, kid(%1%), face(%2%), sibling(%3%)",
            kid % boost::uint16_t(f) % sib);

        children_[kid].set_sibling_push(f, sib, client_from_this());
    }
} // }}}

void octree_server::tie_child_sibling(
    child_index target_kid
  , face target_f
  , octree_client const& target_sib
    )
{ // {{{
    OCTOPUS_ASSERT_FMT_MSG(
        invalid_face > target_f,
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
            siblings_[source_f] != hpx::invalid_id,
            "source face is a boundary, target_kid(%1%), target_face(%2%), "
            "source_kid(%3%), source_face(%4%), target_sibling(%5%)",
            target_kid % boost::uint16_t(target_f) %
            source_kid % boost::uint16_t(source_f) %
            target_sib);

        // Check if we're at a boundary.
        if (children_[target_kid] != hpx::invalid_id)
        {
            children_[target_kid].set_sibling_push
                (target_f, target_sib, client_from_this());

            source_sib = children_[target_kid];
        }

        else
            source_sib = octree_client(amr_boundary); 
    }

    // We established by assertion earlier that siblings_[source_f] is not NULL.
    // FIXME: If source_sib is NULL and that fact is implicitly known by the
    // caller, then this is non-optimal.
    siblings_[source_f].set_child_sibling_push
        (source_kid, source_f, children_[target_kid]);
} // }}}

///////////////////////////////////////////////////////////////////////////////
// Ghost zone communication

// {{{ Description of the algorithm for send_ghost_zone 
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
// }}}

/// 0.) Send out ghost zone data to our siblings. 
/// 1.) Unlock \a l.
/// 2.) Block until out ghost zones for \a phase are ready.
/// 3.) Relock \a l.
void octree_server::communicate_ghost_zones(
    boost::uint64_t phase
  , mutex_type::scoped_lock& l
    )
{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    OCTOPUS_ASSERT_FMT_MSG(
        phase < ghost_zone_queue_.size(),
        "phase (%1%) is greater than the ghost zone queue length (%2%)",
        phase % ghost_zone_queue_.size());


    boost::array<hpx::future<void>, 6> dependencies;

    for (boost::uint64_t i = 0; i < 6; ++i)
    {
        if (siblings_[i].real()
        {
            // 0.) Send out ghost zone data to our siblings. 
            // NOTE: send_ghost_zone_locked is somewhat compute intensive.
            parent.receive_ghost_zone_push(step_, phase, face(i),
                send_ghost_zone_locked(l));

            dependencies.emplace_back(
                ghost_zone_queue_[phase][i].then(
                    boost::bind(&octree_server::add_ghost_zone,
                        this, face(i), _1))); 
        }

        // A boundary, so we pull.
        else
            dependencies.emplace_back(
                siblings_[i].send_ghost_zone_async(face(i)).then(
                    boost::bind(&octree_server::add_ghost_zone,
                        this, face(i), _1))); 
    }

    {
        // 1.) Unlock l.
        hpx::util::unlock_the_lock<mutex_type::scoped_lock> ul(l);

        // 2.) Block until out ghost zones for phase are ready.
        hpx::wait_all(dependencies);

    } // 3.) Relock l.
}

void octree_server::add_ghost_zone(
    face f ///< Bound parameter.
  , hpx::future<vector3d<std::vector<double> > > zone_f
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT(zone_f.is_ready());

    vector3d<std::vector<double> > zone(zone_f.move());

    // The index of the futures in the vector is the face.
    switch (f)
    {
        ///////////////////////////////////////////////////////////////////////
        // X-axis.
        /// for i in [0, BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [BW, GNX - BW)
        ///             U(i, j, k) = sibling[XL].U(GNX - 2 * BW + i, j, k) 
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
            OCTOPUS_ASSERT(zone.y_length() == gnx - 2 * bw); // [BW, GNX - BW)
            OCTOPUS_ASSERT(zone.z_length() == bw);           // [0, BW)

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
            OCTOPUS_ASSERT(zone.y_length() == gnx - 2 * bw); // [BW, GNX - BW)
            OCTOPUS_ASSERT(zone.z_length() == bw);           // [GNX - BW, GNX)

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

void octree_server::receive_ghost_zone(
    boost::uint64_t step ///< For debugging purposes.
  , boost::uint64_t phase 
  , face f ///< Relative to caller.
  , vector3d<std::vector<double> > const& zone
    )
{ // {{{
    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT_MSG(step_ == step,
        "cross-timestep communication occurred, octree is ill-formed");

    OCTOPUS_ASSERT_FMT_MSG(
        phase < ghost_zone_queue_.size(),
        "phase (%1%) is greater than the ghost zone queue length (%2%)",
        phase % ghost_zone_queue_.size());

    OCTOPUS_ASSERT(f != invalid_face);

    // NOTE (wash): boost::move should be safe here, zone is a temporary, even
    // if we're local to the caller. Plus, ATM set_value requires the value to
    // be moved to it.
    ghost_zone_queue_[phase][invert(f)].post(zone);
} // }}} 

vector3d<std::vector<double> > octree_server::send_ghost_zone(
    face f
    )
{ // {{{
    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    send_ghost_zone_locked(f, l);
} // }}}

// Who ya gonna call? Ghostbusters!
vector3d<std::vector<double> > octree_server::send_ghost_zone_locked(
    face f ///< Our direction, relative to the caller.
  , mutex_type::scoped_lock& l
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

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

                        zone(ii, jj, kk) = U_(-gnx + 2 * bw + i, j, k);
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

                        zone(ii, jj, kk) = U_(i, -gnx + 2 * bw + j, k);
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
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [GNX - BW, GNX) */ bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = gnx - bw; k < gnx; ++k)
                    {
                        // Adjusted indices. 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k - (gnx - bw);

                        zone(ii, jj, kk) = U_(i, j, -gnx + 2 * bw + k);
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

// FIXME: Range checking.
boost::array<boost::uint64_t, 3> map_location(
    face f ///< Our direction, relative to the caller.
  , boost::uint64_t i
  , boost::uint64_t j
  , boost::uint64_t k
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;
    bool const reflect_on_z = config().reflect_on_z;

    boost::array<boost::uint64_t, 3> v;

    v[0] = i;
    v[1] = j;
    v[2] = k;

    switch (f)
    {
        case XU:
            v[0] = bw;
            break;
        case XL:
            v[0] = gnx - bw - 1;
            break;
    
        case YU:
            v[1] = bw;
            break;
        case YL:
            v[1] = gnx - bw - 1;
            break;
    
        case ZU:
            v[2] = (reflect_on_z ? gnx - k - 1 : bw);
            break;
        case ZL:
            v[2] = (reflect_on_z ? gnx - k - 1 : gnx - bw - 1);
            break;

        default:
            OCTOPUS_ASSERT(false);
            break;
    }

    return v;
} // }}}

vector3d<std::vector<double> > octree_server::send_mapped_ghost_zone(
    face f ///< Our direction, relative to the caller.
    )
{ // {{{ 
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

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
                        boost::array<boost::uint64_t, 3> v =
                            map_location(f, gnx - 2 * bw + i, j, k);

                        // Adjusted indices (for output ghost zone). 
                        boost::uint64_t const ii = i;
                        boost::uint64_t const jj = j - bw;
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(v);

                        science().enforce_outflow
                            (f, x_face_coords(v[0] + 1, v[1], v[2]));
                    }

            return zone;
        } 

        /// for i in [GNX - BW, GNX)
        ///     for j in [BW, GNX - BW)
        ///         for k in [BW, GNX - BW)
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
                        boost::array<boost::uint64_t, 3> v =
                            map_location(f, -gnx + 2 * bw + i, j, k);

                        // Adjusted indices (for output ghost zone). 
                        boost::uint64_t const ii = i - (gnx - bw);
                        boost::uint64_t const jj = j - bw;
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(v);

                        science().enforce_outflow
                            (f, x_face_coords(v[0], v[1], v[2]));
                    }

            return zone;
        }

        ///////////////////////////////////////////////////////////////////////
        // Y-axis.
        /// for i in [BW, GNX - BW)
        ///     for j in [0, BW)
        ///         for k in [BW, GNX - BW)
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
                        boost::array<boost::uint64_t, 3> v =
                            map_location(f, i, gnx - 2 * bw + j, k);

                        // Adjusted indices (for output ghost zone). 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j;
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(v);

                        science().enforce_outflow
                            (f, y_face_coords(v[0], v[1] + 1, v[2]));
                    }

            return zone;
        } 

        /// for i in [BW, GNX - BW)
        ///     for j in [GNX - BW, GNX)
        ///         for k in [BW, GNX - BW)
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
                        boost::array<boost::uint64_t, 3> v =
                            map_location(f, i, -gnx + 2 * bw + j, k);

                        // Adjusted indices (for output ghost zone). 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - (gnx - bw);
                        boost::uint64_t const kk = k - bw; 

                        zone(ii, jj, kk) = U_(v);

                        science().enforce_outflow
                            (f, y_face_coords(v[0], v[1], v[2]));
                    }

            return zone;
        }

        ///////////////////////////////////////////////////////////////////////
        // Z-axis.
        /// for i in [BW, GNX - BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [0, BW)
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
                        boost::array<boost::uint64_t, 3> v =
                            map_location(f, i, j, gnx - 2 * bw + k);

                        // Adjusted indices (for output ghost zone). 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k;

                        zone(ii, jj, kk) = U_(v);

                        if (config().reflect_on_z)
                            science().reflect_z(zone(ii, jj, kk));
                        else
                            science().enforce_outflow
                                (f, z_face_coords(v[0], v[1], v[2] + 1));
                    }

            return zone;
        } 

        /// for i in [BW, GNX - BW)
        ///     for j in [BW, GNX - BW)
        ///         for k in [GNX - BW, GNX)
        case ZU:
        {
            vector3d<std::vector<double> > zone
                (
                /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [BW, GNX - BW) */  gnx - 2 * bw
              , /* [GNX - BW, GNX) */ bw
                );

            for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
                for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
                    for (boost::uint64_t k = gnx - bw; k < gnx; ++k)
                    {
                        boost::array<boost::uint64_t, 3> v =  
                            map_location(f, i, j, -gnx + 2 * bw + k);

                        // Adjusted indices (for output ghost zone). 
                        boost::uint64_t const ii = i - bw;
                        boost::uint64_t const jj = j - bw; 
                        boost::uint64_t const kk = k - (gnx - bw);

                        zone(ii, jj, kk) = U_(v);

                        if (config().reflect_on_z)
                            science().reflect_z(zone(ii, jj, kk));
                        else
                            science().enforce_outflow
                                (f, z_face_coords(v[0], v[1], v[2]));
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

///////////////////////////////////////////////////////////////////////////////
// Child -> parent injection of state.

/// 0.) Unlock \a l.
/// 1.) Blocks until the child -> parent injection that is at \a phase in 
///     the queue is ready.
/// 2.) Relock \a l.
/// 3.) Send a child -> parent injection up to our parent. 
void octree_server::child_to_parent_injection(
    boost::uint64_t phase
  , mutex_type::scoped_lock& l
    )
{ // {{{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
    
    OCTOPUS_ASSERT_FMT_MSG(
        phase < children_state_queue_.size(),
        "phase (%1%) is greater than the children state queue length (%2%)",
        phase % children_state_queue_.size());

    std::vector<hpx::future<void> > dependencies;
    dependencies.reserve(8);

    bool has_children = false;

    for (boost::uint64_t i = 0; i < 8; ++i)
        if (children_[i] != hpx::invalid_id())
        {
            has_children = true;
            break;
        }

    if (has_children)
    {
        OCTOPUS_ASSERT(config().max_refinement_level != level_)

        for (boost::uint64_t i = 0; i < 8; ++i)
            if (children_[i] != hpx::invalid_id())
                dependencies.emplace_back(
                    children_state_queue_[phase][i].then(
                        boost::bind(&octree_server::add_child_state,
                            this, child_index(i), _1))); 
    }

    {
        // 0.) Unlock l.
        hpx::util::unlock_the_lock<mutex_type::scoped_lock> ul(l);

        // 1.) Blocks until the child -> parent injection that is at phase in 
        //     the queue is ready.
        hpx::wait_all(dependencies);

    } // 2.) Relock l.

    if (parent != hpx::invalid_id)
    {
        OCTOPUS_ASSERT(level_ != 0);

        // 3.) Send a child -> parent injection up to our parent. 
        // NOTE: send_child_state_locked is somewhat compute intensive.
        parent.receive_child_state_push(step_, phase, get_child_index_locked(l),
            send_child_state_locked(l));
    }
} // }}}

void octree_server::add_child_state(
    child_index idx ///< Bound parameter.
  , hpx::future<vector3d<std::vector<double> > > state_f
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    // Make sure that we are initialized.
    initialized_.wait();

    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT(state_f.is_ready());

    vector3d<std::vector<double> > state(state_f.move());

    OCTOPUS_ASSERT(state.x_length() == ((gnx - 2 * bw) / 2));
    OCTOPUS_ASSERT(state.y_length() == ((gnx - 2 * bw) / 2));
    OCTOPUS_ASSERT(state.z_length() == ((gnx - 2 * bw) / 2));

    // NOTE; i += 2 comes from the original code.
    for (boost::uint64_t i = bw; i < (gnx - bw); i += 2)
        for (boost::uint64_t j = bw; j < (gnx - bw); j += 1)
            for (boost::uint64_t k = bw; k < (gnx - bw); k += 1)
            {
                // Adjusted indices (for destination).
                boost::uint64_t const id
                    = (bw + i) / 2 + idx.x() * ((gnx / 2) - bw);
                boost::uint64_t const jd
                    = (bw + j) / 2 + idx.y() * ((gnx / 2) - bw);
                boost::uint64_t const kd
                    = (bw + k) / 2 + idx.z() * ((gnx / 2) - bw);

                // Adjusted indices (for source).
                boost::uint64_t is = i - bw; 
                boost::uint64_t js = j - bw; 
                boost::uint64_t ks = k - bw; 

                U_(id, jd, kd) = state(is, js, ks); 
            }
} // }}}

void octree_server::receive_child_state(
    boost::uint64_t step ///< For debugging purposes.
  , boost::uint64_t phase 
  , child_index idx 
  , vector3d<std::vector<double> > const& state
    )
{ // {{{
    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT_MSG(step_ == step,
        "cross-timestep communication occurred, octree is ill-formed");

    OCTOPUS_ASSERT_FMT_MSG(
        phase < children_state_queue_.size(),
        "phase (%1%) is greater than the children state queue length (%2%)",
        phase % children_state_queue_.size());

    OCTOPUS_ASSERT(boost::uint64_t(idx) < 8);

    // NOTE (wash): boost::move should be safe here, zone is a temporary, even
    // if we're local to the caller. Plus, ATM set_value requires the value to
    // be moved to it.
    children_state_queue_[phase][idx].post(state);
} // }}}

vector3d<std::vector<double> > send_child_state_locked(
    mutex_type::scoped_lock &l
    )
{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    vector3d<std::vector<double> > state
        (
        (gnx - 2 * bw) / 2
      , (gnx - 2 * bw) / 2
      , (gnx - 2 * bw) / 2
        ); 

    child_index c = get_child_index_locked(l);

    // NOTE; i += 2 comes from the original code.
    for (boost::uint64_t i = bw; i < (gnx - bw); i += 2)
        for (boost::uint64_t j = bw; j < (gnx - bw); j += 1)
            for (boost::uint64_t k = bw; k < (gnx - bw); k += 1)
            {
                // Adjusted indices (for output state).
                boost::uint64_t const ii = ((i + bw) / 2) - bw;
                boost::uint64_t const jj = ((j + bw) / 2) - bw; 
                boost::uint64_t const kk = ((k + bw) / 2) - bw; 

                state(ii, jj, kk) += U_(i + 0, j + 0, k + 0);
                state(ii, jj, kk) += U_(i + 1, j + 0, k + 0);
                state(ii, jj, kk) += U_(i + 0, j + 1, k + 0);
                state(ii, jj, kk) += U_(i + 1, j + 1, k + 0);
                state(ii, jj, kk) += U_(i + 0, j + 0, k + 1);
                state(ii, jj, kk) += U_(i + 1, j + 0, k + 1);
                state(ii, jj, kk) += U_(i + 0, j + 1, k + 1);
                state(ii, jj, kk) += U_(i + 1, j + 1, k + 1);
                state(ii, jj, kk) *= 0.125;
            }

    return state; 
}

///////////////////////////////////////////////////////////////////////////////
// Tree traversal.
void octree_server::apply(
    hpx::util::function<void(octree_server&)> const& f
    )
{ // {{{
    std::vector<hpx::future<void> > recursion_is_parallelism;

    // Make sure that we are initialized.
    initialized_.wait();
    
    recursion_is_parallelism.reserve(8);
    
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
            recursion_is_parallelism.push_back(children_[i].apply_async(f)); 
    
    // Invoke the kernel on ourselves ...
    f(*this);

    // ... and block while our children compute.
    hpx::wait(recursion_is_parallelism); 
} // }}}

void octree_server::apply_leaf(
    hpx::util::function<void(octree_server&)> const& f
    )
{ // {{{
    std::vector<hpx::future<void> > recursion_is_parallelism;

    // Make sure that we are initialized.
    initialized_.wait();
    
    f(*this);
} // }}}

///////////////////////////////////////////////////////////////////////////////

// The "boring" version. Does one step. 
void octree_server::step(double dt)
{ // {{{
    OCTOPUS_ASSERT_MSG(0 < dt, "invalid timestep size");

    std::vector<hpx::future<void> > recursion_is_parallelism;

    {
        // Make sure that we are initialized.
        initialized_.wait();
    
        mutex_type::scoped_lock l(mtx_);
    
        recursion_is_parallelism.reserve(8);
    
        // Start recursively executing the kernel function on our children.
        for (boost::uint64_t i = 0; i < 8; ++i)
            if (hpx::invalid_id != children_[i])
                recursion_is_parallelism.push_back(children_[i].step_async(dt)); 
    
        // Kernel.
        step_kernel(dt, l);
    }

    // Block while our children compute.
    hpx::wait(recursion_is_parallelism); 
} // }}}

// Recursion is parallel iteration.
// IMPLEMENT (properly).
void octree_server::step_to_time(double dt, double until)
{ // {{{
    OCTOPUS_ASSERT_MSG(0 < dt, "invalid timestep size");

/*
    // Make sure that we are initialized.
    initialized_.wait();

    {
        mutex_type::scoped_lock l(mtx_);

        // Start recursively executing the kernel on our children.
        for (boost::uint64_t i = 0; i < 8; ++i)
            if (hpx::invalid_id != children_[i])
                children_[i].step_to_time_push(dt, until); 

        ///////////////////////////////////////////////////////////////////////
        // Kernel.
        step_kernel(dt, l);
    }

    // Are we done?
    if ((level_ == 0) && (until <= time_))
    {
        // If not, keep going. 
        double next_dt = science().next_timestep(*this, step_, dt, until);

        OCTOPUS_ASSERT(future_self != hpx::invalid_id);

        // More recursion!
        future_self.step_to_time_push(next_dt, until);
    }
*/
} // }}}

void octree_server::step_kernel(
    double dt
  , mutex_type::scoped_lock& l
    )
{ // {{{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    U0_ = U_;
    FO0_ = FO_;

    // NOTE (to self): I have no good place to put this, so I'm putting it
    // here: we do TVD RK3 (google is your friend).
    switch (config().runge_kutta_order)
    {
        case 1:
        {
            sub_step_kernel(0, dt, 1.0, l);
            break;
        }

        case 2:
        {
            sub_step_kernel(0, dt, 1.0, l);
            sub_step_kernel(1, dt, 0.5, l);
            break; 
        }

        case 3:
        {
            sub_step_kernel(0, dt, 1.0, l);
            sub_step_kernel(1, dt, 0.25, l);
            sub_step_kernel(2, dt, 2.0 / 3.0, l);
            break; 
        }

        default:
        {
            OCTOPUS_ASSERT_FMT_MSG(false,
                "runge-kutta order (%1%) is unsupported or invalid",
                config().runge_kutta_order);
        }
    };

    // Talks, unlocks, blocks, and relocks.
    // REVIEW (wash): Is this /really/ necessary? We do ghost zones at the
    // beginning of the next timestep already.
    // NOTE (wash): See math in the header above the declaration of
    // ghost_zone_queue for an explaination of the phase used here.
    communicate_ghost_zones(config().runge_kutta_order, l);

    ++step_;
    time_ += dt;
} // }}}

// Two communication phases.
void octree_server::sub_step_kernel(
    boost::uint64_t phase
  , double dt
  , double beta
  , mutex_type::scoped_lock &l 
    )
{ // {{{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    // Talks, unlocks, blocks, and relocks.
    communicate_ghost_zones(phase, l);

    // Semi-trivial.
    prepare_differentials_kernel(l);

    // Operations parallelizes by axis.
    compute_flux_kernel(l);
    adjust_flux_kernel(l);

    // Critical section.
    sum_differentials_kernel(l);
    add_differentials_kernel(dt, beta, l);

    // Unlocks, blocks, relocks and talks. 
    child_to_parent_injection(phase, l);
} // }}}

void octree_server::add_differentials_kernel(
    double dt
  , double beta
  , mutex_type::scoped_lock &l 
    )
{ // {{{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    using namespace octopus::operators;

    for (boost::uint64_t i = bw; i < gnx - bw; ++i)
        for (boost::uint64_t j = bw; j < gnx - bw; ++j)
            for (boost::uint64_t k = bw; k < gnx - bw; ++k)
            {
                boost::array<double, 3> coords = center_coords(i, j, k);

                D_(i, j, k) += science().source(U_(i, j, k), coords);

                // Here, you can see the temporal dependency.
                U_(i, j, k) = (U_(i, j, k) + D_(i, j, k) * dt) * beta
                            + U0_(i, j, k) * (1.0 - beta); 

                science().floor(U_(i, j, k), coords);
            }

    FO_ = (FO_ + DFO_ * dt) * beta + FO0_ * (1.0 - beta);
} // }}}

void octree_server::prepare_differentials_kernel(
    mutex_type::scoped_lock &l
    ) 
{ // {{{
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const gnx = config().grid_node_length;

    OCTOPUS_ASSERT(DFO_.size() == ss);

    for (boost::uint64_t i = 0; i < ss; ++i)
        DFO_[i] = 0.0;

    OCTOPUS_ASSERT(D_.size() == (gnx * gnx * gnx));

    for (boost::uint64_t i = 0; i < gnx; ++i)
        for (boost::uint64_t j = 0; j < gnx; ++j)
            for (boost::uint64_t k = 0; k < gnx; ++k)
                for (boost::uint64_t l = 0; l < ss; ++l)
                {
                    OCTOPUS_ASSERT(D_(i, j, k).size() == ss);
                    D_(i, j, k)[l] = 0.0;
                }
} // }}}

void octree_server::compute_flux_kernel(
    mutex_type::scoped_lock &l 
    )
{ // {{{ 
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    ////////////////////////////////////////////////////////////////////////////    
    // Compute our own local fluxes locally in parallel. 

    //hpx::util::unlock_the_lock<mutex_type::scoped_lock> ul(l);

    // Do two in other threads.
    boost::array<hpx::future<void>, 2> xyz =
    { {
        hpx::async(boost::bind
            (&octree_server::compute_x_flux_kernel, this, boost::ref(l)))
      , hpx::async(boost::bind
            (&octree_server::compute_y_flux_kernel, this, boost::ref(l)))
    } };

    // And do one here.
    compute_z_flux_kernel(l);

    // Wait for the local x and y fluxes to be computed.
    hpx::wait(xyz[0], xyz[1]);
} // }}}

void octree_server::compute_x_flux_kernel(
    mutex_type::scoped_lock &l 
    )
{ // {{{ 
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    std::vector<std::vector<double> > q0(gnx, std::vector<double>(ss));
    std::vector<std::vector<double> > ql(gnx, std::vector<double>(ss));
    std::vector<std::vector<double> > qr(gnx, std::vector<double>(ss));

    for (boost::uint64_t k = bw; k < (gnx - bw); ++k)
        for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
        {
            for (boost::uint64_t i = 0; i < gnx; ++i)
            {
                q0[i] = U_(i, j, k);
    
                boost::array<double, 3> coords = center_coords(i, j, k);
    
                science().conserved_to_primitive(q0[i], coords);
            }
    
            science().reconstruct(q0, ql, qr);
    
            for (boost::uint64_t i = bw; i < gnx - bw + 1; ++i)
            {
                boost::array<double, 3> coords = x_face_coords(i, j, k);
    
                science().primitive_to_conserved(ql[i], coords);
                science().primitive_to_conserved(qr[i], coords);
    
                double const a =
                    (std::max)(science().max_eigenvalue(x_axis, ql[i], coords)
                             , science().max_eigenvalue(x_axis, qr[i], coords));
    
                std::vector<double>
                    ql_flux = science().flux(x_axis, ql[i], coords),
                    qr_flux = science().flux(x_axis, qr[i], coords);
    
                using namespace octopus::operators;
     
                FX_(i, j, k) = ((ql_flux + qr_flux) - (qr[i] - ql[i]) * a)
                             * 0.5;
            }
        }
} // }}}

void octree_server::compute_y_flux_kernel(
    mutex_type::scoped_lock &l 
    )
{ // {{{ 
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    std::vector<std::vector<double> > q0(gnx, std::vector<double>(ss));
    std::vector<std::vector<double> > ql(gnx, std::vector<double>(ss));
    std::vector<std::vector<double> > qr(gnx, std::vector<double>(ss));

    for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
        for (boost::uint64_t k = bw; k < (gnx - bw); ++k)
        {
            for (boost::uint64_t j = 0; j < gnx; ++j)
            {
                q0[j] = U_(i, j, k);
    
                boost::array<double, 3> coords = center_coords(i, j, k);
    
                science().conserved_to_primitive(q0[j], coords);
            }
    
            science().reconstruct(q0, ql, qr);
    
            for (boost::uint64_t j = bw; j < gnx - bw + 1; ++j)
            {
                boost::array<double, 3> coords = y_face_coords(i, j, k);
    
                science().primitive_to_conserved(ql[j], coords);
                science().primitive_to_conserved(qr[j], coords);
    
                double const a =
                    (std::max)(science().max_eigenvalue(y_axis, ql[j], coords)
                             , science().max_eigenvalue(y_axis, qr[j], coords));
    
                std::vector<double>
                    ql_flux = science().flux(y_axis, ql[j], coords)
                  , qr_flux = science().flux(y_axis, qr[j], coords);
     
                using namespace octopus::operators;
    
                FY_(i, j, k) = ((ql_flux + qr_flux) - (qr[j] - ql[j]) * a)
                             * 0.5;
            }
        }
} // }}}

void octree_server::compute_z_flux_kernel(
    mutex_type::scoped_lock &l 
    )
{ // {{{ 
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const ss = science().state_size;
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    indexer2d<2> const indexer(bw, gnx - bw - 1, bw, gnx - bw - 1);

    std::vector<std::vector<double> > q0(gnx, std::vector<double>(ss));
    std::vector<std::vector<double> > ql(gnx, std::vector<double>(ss));
    std::vector<std::vector<double> > qr(gnx, std::vector<double>(ss));

    for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
        for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
        {
            for (boost::uint64_t k = 0; k < gnx; ++k)
            {
                q0[k] = U_(i, j, k);
    
                boost::array<double, 3> coords = center_coords(i, j, k);

                science().conserved_to_primitive(q0[k], coords);
            }
    
            science().reconstruct(q0, ql, qr);
    
            for (boost::uint64_t k = bw; k < gnx - bw + 1; ++k)
            {
                boost::array<double, 3> coords = z_face_coords(i, j, k);
    
                science().primitive_to_conserved(ql[k], coords);
                science().primitive_to_conserved(qr[k], coords);
    
                double const a =
                    (std::max)(science().max_eigenvalue(z_axis, ql[k], coords)
                             , science().max_eigenvalue(z_axis, qr[k], coords));
    
                std::vector<double>
                    ql_flux = science().flux(z_axis, ql[k], coords)
                  , qr_flux = science().flux(z_axis, qr[k], coords)
                    ;
     
                using namespace octopus::operators;
    
                FZ_(i, j, k) = ((ql_flux + qr_flux) - (qr[k] - ql[k]) * a) * 0.5;
            }
        }
} // }}}

void octree_server::adjust_flux_kernel(
    mutex_type::scoped_lock& l
    )
{ // {{{ IMPLEMENT
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

} // }}}

void octree_server::sum_differentials_kernel(
    mutex_type::scoped_lock& l
    )
{ // {{{ 
    OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");

    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    double const dxinv = 1.0 / dx_;

    ///////////////////////////////////////////////////////////////////////////
    // Kernel.
    using namespace octopus::operators;

    // NOTE: This is probably too tight a loop to parallelize with HPX, but
    // could be vectorized. 
    for (boost::uint64_t i = bw; i < gnx - bw; ++i)
        for (boost::uint64_t j = bw; j < gnx - bw; ++j)
        {
            for (boost::uint64_t k = bw; k < gnx - bw; ++k)
            {
                D_(i, j, k) -= (FX_(i + 1, j, k) - FX_(i, j, k)) * dxinv;
                D_(i, j, k) -= (FY_(i, j + 1, k) - FY_(i, j, k)) * dxinv;
                D_(i, j, k) -= (FZ_(i, j, k + 1) - FZ_(i, j, k)) * dxinv;
            }

            // Only do this for the root node.
            // REVIEW: Talk to Dominic about how this works, this is ported
            // directly from the original code.
            if (level_ == 0)
            {
                // i = y-axis, j = z-axis 
                DFO_ += (FX_(gnx - bw, i, j) - FX_(bw, i, j)) * dx_ * dx_;

                // i = x-axis, j = z-axis 
                DFO_ += (FY_(i, gnx - bw, j) - FY_(i, bw, j)) * dx_ * dx_;
        
                // i = x-axis, j = y-axis 
                if (config().reflect_on_z)
                    DFO_ += (FZ_(i, j, gnx - bw)) * dx_ * dx_;
                else
                    DFO_ += (FZ_(i, j, gnx - bw) - FZ_(i, j, bw)) * dx_ * dx_;
            }
        }
} // }}}

void octree_server::copy_and_regrid()
{ // {{{ IMPLEMENT

} // }}}

void octree_server::output()
{ // {{{ 
    std::vector<hpx::future<void> > recursion_is_parallelism;

    // Make sure that we are initialized.
    initialized_.wait();
    
    mutex_type::scoped_lock l(mtx_);
    
    recursion_is_parallelism.reserve(8);
    
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
            recursion_is_parallelism.push_back(
                hpx::async<output_action>(children_[i].get_gid()));

    science().output(*this); 

    hpx::wait(recursion_is_parallelism); 
} // }}}

}

