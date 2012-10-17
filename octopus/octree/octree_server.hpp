////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0)
#define OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

#include <hpx/runtime/components/server/managed_component_base.hpp>
#include <hpx/lcos/local/mutex.hpp>
#include <hpx/lcos/local/event.hpp>

#include <octopus/octree/octree_init_data.hpp>
#include <octopus/octree/octree_client.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT octree_server
  : hpx::components::managed_component_base<octree_server>
{
  private:
    typedef hpx::lcos::local::mutex mutex_type;

    hpx::lcos::local::event initialized_;
    mutable mutex_type mtx_; 
    boost::uint8_t siblings_set_;
    bool received_state_;

    ///////////////////////////////////////////////////////////////////////////
    // From OctNode
    octree_client parent_; 
    boost::array<octree_client, 8> children_;
    boost::array<octree_client, 6> siblings_; 
    boost::uint64_t level_;
    boost::array<boost::int64_t, 3> location_; ///< NOTE: Not sure if this needs
                                               /// to be signed.

    ///////////////////////////////////////////////////////////////////////////
    // From Grid/GridNode
    double dx_;   ///< The spatial size of the node (w/o ghost zones).
                  ///  NOTE: Confirmation needed from Dominic. 
    double time_; ///< The current (physics?) time.
                  ///  NOTE: Confirmation needed from Dominic.

    boost::array<boost::int64_t, 3> offset_; ///< NOTE: Not sure if this needs
                                             ///  be signed. 

    boost::array<double, 3> origin_; ///< The origin of the cartesian grid.
                                     ///  NOTE: Confirmation needed from
                                     ///  from Dominic.

    boost::uint64_t step_;

    // REVIEW: Consider compile-time maximum sizes for the state vector, to
    // optimize allocations.
    // REVIEW: Are the "corners" of our cube are needed? I think there may be
    // eight chunks of unused space that we don't ever use.
    vector3d<std::vector<double> > U_; ///< 3d array of state vectors, includes
                                       ///  ghost zones. Size of the state
                                       ///  vectors comes from the science
                                       ///  table.

    vector3d<std::vector<double> > U0_; ///< 3d array of state vectors, includes
                                        ///  ghost zones. Used for summing the
                                        ///  differential. 

    ///////////////////////////////////////////////////////////////////////////
    // TODO: Migration.
#if 0
    friend class boost::serialization::access;

    // FIXME: Should not be able to move if not initialized.
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        // IMPLEMENT
    }
#endif

    // Precondition: mtx_ must be locked.
    child_index get_child_index_locked(mutex_type::scoped_lock& l) const
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(0 != level_, "root octree_server has no parent");
        child_index idx(location_[0] % 2, location_[1] % 2, location_[2] % 2);
        return idx; 
    }

    // Precondition: mtx_ must be locked.
    void sibling_set_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(siblings_set_ < 6, "double initialization");

        for (std::size_t i = 0; i < 6; ++i)
            OCTOPUS_ASSERT(siblings_[i] != hpx::invalid_id);

        if ((++siblings_set_ == 6) && received_state_)
            initialized_.set(); 
    }  

    // Precondition: mtx_ must be locked.
    void state_received_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(received_state_ == true, "double initialization");
        received_state_ = true;
        if (siblings_set_ == 6)
            initialized_.set(); 
    }  

    child_index get_child_index() const
    {
        mutex_type::scoped_lock l(mtx_);
        return get_child_index_locked(l); 
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set the state of this node based on the state of its parent. 
    /// 
    /// Remote Operations:   No.
    /// Concurrency Control: Locks mtx_.
    /// Synchrony Gurantee:  Synchronous. 
    void inject_state_from_parent(
        vector3d<std::vector<double> > const& pU
        );

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Get a reference to this node that is safe to pass to our
    ///        children. Said reference must be uncounted to prevent reference
    ///        cycles.
    /// 
    /// Remote Operations:   No.
    /// Concurrency Control: No.
    /// Synchrony Gurantee:  Synchronous. 
    hpx::id_type reference_from_this()
    {
        // We shouldn't need to lock here, I believe.
        hpx::id_type gid = get_gid();
        OCTOPUS_ASSERT_MSG(
            gid.get_management_type() == hpx::id_type::unmanaged,
            "get_gid() should return an unmanaged GID");
        return gid;
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Get a client referencing this node that is safe to pass to our
    ///        children. The reference the child holds must be uncounted to
    ///        prevent reference cycles.
    /// 
    /// Remote Operations:   No.
    /// Concurrency Control: No.
    /// Synchrony Gurantee:  Synchronous. 
    octree_client client_from_this()
    {
        // We shouldn't need to lock here, I believe.
        hpx::id_type gid = get_gid();
        OCTOPUS_ASSERT_MSG(
            gid.get_management_type() == hpx::id_type::unmanaged,
            "get_gid() should return an unmanaged GID");
        return octree_client(gid);;
    }

  public:
    octree_server()
    {
        OCTOPUS_ASSERT_MSG(false, "octree_server can't be default constructed");
    } 

    // FIXME: Create physical bounds?
    /// \brief Construct a root node 
    octree_server(
        octree_init_data const& init
        )
      : siblings_set_(6)
      , received_state_(false)
      , parent_(init.parent)
      , level_(init.level)
      , location_(init.location)
      , dx_(init.dx)
      , time_(init.time)
      , offset_(init.offset)
      , origin_(init.origin)
      , step_(0)
    {
        OCTOPUS_TEST_IN_PLACE(parent_ == hpx::invalid_id);
    }

    // FIXME: Non-optimal, inject_state_from_parent should be called with
    // fire-and-forget semantics. However, the lifetime of parent_U then
    // becomes an issue. Because of this, U_ may need to be made into a
    // shared_ptr.
    /// \brief Construct a child node.
    octree_server(
        octree_init_data const& init
      , vector3d<std::vector<double> > const& parent_U
        )
      : siblings_set_(0)
      , received_state_(false)
      , parent_(init.parent)
      , level_(init.level)
      , location_(init.location)
      , dx_(init.dx)
      , time_(init.time)
      , offset_(init.offset)
      , origin_(init.origin)
      , step_(0)
    {
        // Make sure our parent reference is not reference counted.
        OCTOPUS_ASSERT_MSG(
            init.parent.get_management_type() == hpx::id_type::unmanaged,
            "reference cycle detected in child");
 
        inject_state_from_parent(parent_U);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Create the \a kid child for this node.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantee:  Fire-and-Forget 
    void create_child(
        child_index kid
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                create_child,
                                create_child_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set \a target_sib as the \a target_f sibling of this node.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Locks mtx_.
    /// Synchrony Gurantee:  Fire-and-Forget 
    void set_sibling(
        face f
      , octree_client const& sib
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                set_sibling,
                                set_sibling_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set \a target_sib as the \a target_f sibling of this node.
    ///        Additionally, set this node as the invert(target_f) sibling of
    ///        \a target_sib.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Locks mtx_.
    /// Synchrony Gurantee:  Fire-and-Forget. 
    void tie_sibling(
        face target_f
      , octree_client const& target_sib
      , octree_client const& target_sib_parent
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                tie_sibling,
                                tie_sibling_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set \a target_sib as the \a target_f sibling of this node's
    ///        \a target_kid child.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantee:  Fire-and-Forget. 
    void set_child_sibling(
        child_index kid
      , face f
      , octree_client const& sib
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                set_child_sibling,
                                set_child_sibling_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set \a target_sib as the \a target_f sibling of this node's
    ///        \a target_kid child. Additionally, set this node's \a target_kid
    ///        child as the invert(target_f) sibling of \a target_sib.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantee:  Fire-and-Forget.
    void tie_child_sibling(
        child_index target_kid
      , face target_f
      , octree_client const& target_sib
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                tie_child_sibling,
                                tie_child_sibling_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set \a target_sib as the \a target_f sibling of this node's
    ///        \a target_kid child. Additionally, set this node's \a target_kid
    ///        child as the invert(target_f) sibling of \a target_sib.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantee:  Fire-and-Forget.
    boost::array<octree_client, 6> get_siblings();

    // REVIEW: Should this be a direct action?
    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                get_siblings,
                                get_siblings_action);

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT
    void inject_state_from_children();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                inject_state_from_children,
                                inject_state_from_children_action);


    ///////////////////////////////////////////////////////////////////////////
    // FIXME: Push don't pull.
    // NOTE: enforce_boundaries in the original code.
    /// \brief Requests ghost zone data from all siblings. 
    void receive_ghost_zones();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                receive_ghost_zones,
                                receive_ghost_zones_action);

    // FIXME: Push don't pull.
    /// \brief Produces ghost zone data for a sibling.
    vector3d<std::vector<double> > send_ghost_zone(
        face f ///< Our direction, relative to the caller.
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                send_ghost_zone,
                                send_ghost_zone_action);

  private:
    // FIXME: Rvalue reference kung-fo must be applied here
    void integrate_ghost_zone(
        std::size_t i
      , vector3d<std::vector<double> > const& zone
        );

  public:

    ///////////////////////////////////////////////////////////////////////////
    // NOTE: exec_function in the original code.
    void apply(
        hpx::util::function<void(octree_server&)> const& f
      , boost::uint64_t minimum_level
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                apply,
                                apply_action);

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT
    void save_state();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                save_state,
                                save_state_action);  

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT
    void add_differentials(double dt, double beta); 

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                add_differentials,
                                add_differentials_action);  

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT
    void clear_differentials(); 

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                clear_differentials,
                                clear_differentials_action);  

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Evolve the system for a temporal period of \a dt.
    // NOTE: The implementation of this should be in the science table, or
    // maybe the whole thing should live in the driver.
    // NOTE: This function DOES NOT LOCK, do not call concurrently without
    // synchronization. 
    void step(double dt);

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                step,
                                step_action);  

  private:
    // NOTE: The implementation of this should be in the science table, or
    // maybe the whole thing should live in the driver.
    // NOTE: This function DOES NOT LOCK, do not call concurrently without
    // synchronization.
    void sub_step(double dt, double beta);

  public:

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT
    void refine();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                refine,
                                refine_action);  

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT 
    void compute_x_flux();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                compute_x_flux,
                                compute_x_flux_action);  

    // IMPLEMENT
    void compute_y_flux();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                compute_y_flux,
                                compute_y_flux_action);  

    // IMPLEMENT
    void compute_z_flux();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                compute_z_flux,
                                compute_z_flux_action);  

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT 
    void adjust_x_flux();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                adjust_x_flux,
                                adjust_x_flux_action);  

    // IMPLEMENT
    void adjust_y_flux();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                adjust_y_flux,
                                adjust_y_flux_action);  

    // IMPLEMENT
    void adjust_z_flux();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                adjust_z_flux,
                                adjust_z_flux_action);  

    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT 
    void sum_x_differentials();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                sum_x_differentials,
                                sum_x_differentials_action);  

    // IMPLEMENT
    void sum_y_differentials();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                sum_y_differentials,
                                sum_y_differentials_action);  

    // IMPLEMENT
    void sum_z_differentials();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                sum_z_differentials,
                                sum_z_differentials_action);  
};

}

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION_DECLARATION(                                        \
        octopus::octree_server::BOOST_PP_CAT(name, _action),                \
        BOOST_PP_CAT(octopus_octree_server_, BOOST_PP_CAT(name, _action)))  \
    /**/

OCTOPUS_REGISTER_ACTION(create_child);
OCTOPUS_REGISTER_ACTION(set_sibling);
OCTOPUS_REGISTER_ACTION(tie_sibling);
OCTOPUS_REGISTER_ACTION(set_child_sibling);
OCTOPUS_REGISTER_ACTION(tie_child_sibling);
OCTOPUS_REGISTER_ACTION(get_siblings);
OCTOPUS_REGISTER_ACTION(inject_state_from_children);
OCTOPUS_REGISTER_ACTION(send_ghost_zone);
OCTOPUS_REGISTER_ACTION(receive_ghost_zones);
OCTOPUS_REGISTER_ACTION(apply);
OCTOPUS_REGISTER_ACTION(save_state);
OCTOPUS_REGISTER_ACTION(add_differentials);
OCTOPUS_REGISTER_ACTION(clear_differentials);
OCTOPUS_REGISTER_ACTION(step);
OCTOPUS_REGISTER_ACTION(refine);
OCTOPUS_REGISTER_ACTION(compute_x_flux);
OCTOPUS_REGISTER_ACTION(compute_y_flux);
OCTOPUS_REGISTER_ACTION(compute_z_flux);
OCTOPUS_REGISTER_ACTION(adjust_x_flux);
OCTOPUS_REGISTER_ACTION(adjust_y_flux);
OCTOPUS_REGISTER_ACTION(adjust_z_flux);
OCTOPUS_REGISTER_ACTION(sum_x_differentials);
OCTOPUS_REGISTER_ACTION(sum_y_differentials);
OCTOPUS_REGISTER_ACTION(sum_z_differentials);

#undef OCTOPUS_REGISTER_ACTION

#endif // OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

