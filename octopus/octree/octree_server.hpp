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

// This is either really beautiful, really messy, or both.

namespace octopus
{

struct OCTOPUS_EXPORT octree_server; 

}

namespace hpx { namespace traits
{

template <>
struct managed_component_ctor_policy<octopus::octree_server>
{
    typedef hpx::traits::construct_with_back_ptr type;
};

}}

namespace octopus
{

struct OCTOPUS_EXPORT octree_server
  : hpx::components::managed_component_base<octree_server>
{
  private:
    friend struct checkout_state;

    typedef hpx::components::managed_component_base<octree_server> base_type;

    typedef hpx::components::managed_component<octree_server>*
        back_pointer_type;

    typedef hpx::lcos::local::mutex mutex_type;

    ///< This event is triggered when we are first initialized. 
    hpx::lcos::local::event initialized_;

    mutable mutex_type mtx_; 
    boost::uint8_t siblings_set_;
    bool state_received_;

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
    // NOTE: Rename dx_;
    double dx_;   ///< The spatial size of the node (w/o ghost zones). h in the
                  ///  original code. NOTE: Confirmation needed from Dominic.  

    double dx0_;

    double time_; ///< The current (physics?) time.
                  ///  NOTE: Confirmation needed from Dominic.

    boost::array<boost::int64_t, 3> offset_; ///< NOTE: Not sure if this needs
                                             ///  be signed. 

    boost::array<double, 3> origin_; ///< The origin of the cartesian grid.
                                     ///  NOTE: Confirmation needed from
                                     ///  from Dominic.

    // TODO: Rename step_.
    boost::uint64_t step_;

    // REVIEW: Consider compile-time maximum sizes for the state vector, to
    // optimize allocations.
    // REVIEW: Are the "corners" of our cube are needed? I think there may be
    // eight chunks of unused space that we don't ever use.
    vector3d<std::vector<double> > U_; ///< 3d array of state vectors, includes
                                       ///  ghost zones. Size of the state
                                       ///  vectors comes from the science
                                       ///  table.

    vector3d<std::vector<double> > U0_; ///< Data from the previous timestep.

    vector3d<std::vector<double> > FX_; ///< Flux (X-axis).
    vector3d<std::vector<double> > FY_; ///< Flux (Y-axis).
    vector3d<std::vector<double> > FZ_; ///< Flux (Z-axis).

    std::vector<double> FO_; ///< Flow off (stuff that leaves the problem
                             ///  space).

    std::vector<double> FO0_;

    vector3d<std::vector<double> > D_; ///< Flux differential.

    std::vector<double> DFO_; ///< Flow off differential. 

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

    // Preconditions: mtx_ must be locked, siblings_set_ must be less than 6.
    void sibling_set_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(siblings_set_ < 6, "double initialization");

        for (std::size_t i = 0; i < 6; ++i)
            OCTOPUS_ASSERT(siblings_[i] != hpx::invalid_id);

        if ((++siblings_set_ == 6) && state_received_)
            initialized_.set(); 
    }  

    // Preconditions: mtx_ must be locked, state_received_ is false. 
    void state_received_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(state_received_ == false, "double initialization");

        state_received_ = true;

        for (std::size_t i = 0; i < 6; ++i)
            OCTOPUS_ASSERT(siblings_[i] != hpx::invalid_id);

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
    octree_server(back_pointer_type back_ptr)
      : base_type(back_ptr) 
    {
        OCTOPUS_ASSERT_MSG(false, "octree_server can't be default constructed");
    } 

    /// \brief Construct a root node. 
    octree_server(
        back_pointer_type back_ptr
      , octree_init_data const& init
        );

    // FIXME: Non-optimal, inject_state_from_parent should be called with
    // fire-and-forget semantics. However, the lifetime of parent_U then
    // becomes an issue. Because of this, U_ may need to be made into a
    // shared_ptr.
    /// \brief Construct a child node.
    octree_server(
        back_pointer_type back_ptr
      , octree_init_data const& init
      , vector3d<std::vector<double> > const& parent_U
        );

    // FIXME: More descriptive name.
    boost::array<double, 3> xfx(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
    	boost::array<double, 3> x = { { xf(i), yc(j), zc(k) } };
    	return x;
    }

    boost::uint64_t get_level() const
    {
        return level_;
    }

    double get_time() const
    {
        return time_;
    }

    boost::uint64_t get_step() const
    {
        return step_;
    }

    // Only used as an action by output().
    HPX_DEFINE_COMPONENT_CONST_DIRECT_ACTION(octree_server,
                                             get_step,
                                             get_step_action);  

    std::vector<double>& operator()(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        )
    {
        return U_(i, j, k);
    } 

    std::vector<double> const& operator()(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
        return U_(i, j, k);
    } 

    // NOTE: Use with caution.
    mutex_type& get_mutex() const
    {
        return mtx_;
    }

    // FIXME: More descriptive name.
    boost::array<double, 3> xfy(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
    	boost::array<double, 3> x = { { xc(i), yf(j), zc(k) } };
    	return x;
    }
    
    // FIXME: More descriptive name.
    boost::array<double, 3> xfz(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
    	boost::array<double, 3> x = { { xc(i), yc(j), zf(k) } };
    	return x;
    }

    // FIXME: More descriptive name.
    double xc(boost::uint64_t i) const
    {
        return xf(i) + 0.5 * dx_;
    }

    // FIXME: More descriptive name.
    double yc(boost::uint64_t j) const
    {
        return yf(j) + 0.5 * dx_;
    }

    // FIXME: More descriptive name.
    double zc(boost::uint64_t k) const
    {
        return zf(k) + 0.5 * dx_;
    }

    double xf(boost::uint64_t i) const;

    double yf(boost::uint64_t i) const;

    double zf(boost::uint64_t i) const;

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
        face f ///< Caller's direction, relative to us.
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
    boost::array<octree_client, 6> get_siblings() 
    {
        // Make sure that we are initialized.
        initialized_.wait();

        // Is the lock needed?
        mutex_type::scoped_lock l(mtx_);
        return siblings_;
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                get_siblings,
                                get_siblings_action);

    ///////////////////////////////////////////////////////////////////////////
    boost::array<boost::int64_t, 3> get_offset() const
    {
        return offset_;
    }

    HPX_DEFINE_COMPONENT_CONST_DIRECT_ACTION(octree_server,
                                             get_offset,
                                             get_offset_action);

    ///////////////////////////////////////////////////////////////////////////
    // NOTE: enforce_boundaries in the original code.
    // IMPLEMENT: Push don't pull.
    /// \brief Requests ghost zone data from all siblings. 
    void receive_ghost_zones();

  private:
    void receive_ghost_zones_kernel(
        mutex_type::scoped_lock& l
        );

  public:
    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                receive_ghost_zones,
                                receive_ghost_zones_action);

    // IMPLEMENT: Push don't pull.
    /// \brief Produces ghost zone data for a sibling.
    vector3d<std::vector<double> > send_ghost_zone(
        face f ///< Our direction, relative to the caller.
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                send_ghost_zone,
                                send_ghost_zone_action);

    // IMPLEMENT: Push don't pull.
    vector3d<std::vector<double> > send_mapped_ghost_zone(
        face f ///< Our direction, relative to the caller.
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                send_mapped_ghost_zone,
                                send_mapped_ghost_zone_action);

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
      , boost::uint64_t minimum_level = 0
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                apply,
                                apply_action);

    ///////////////////////////////////////////////////////////////////////////
    void step(double dt);

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                step,
                                step_action);  

    /// Perform a timestep of period \a dt. Then, call step_to_time() so long as
    /// \a time_ is less than \a until. 
    void step_to_time(double dt, double until);

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                step_to_time,
                                step_to_time_action);  

  private:
    void step_kernel(
        double dt
      , mutex_type::scoped_lock& l
        );

    void sub_step_kernel(
        double dt
      , double beta
      , mutex_type::scoped_lock& l
        );

    // IMPLEMENT
    void receive_state_from_children_kernel(mutex_type::scoped_lock& l);

    void add_differentials_kernel(
        double dt
      , double beta
      , mutex_type::scoped_lock& l
        ); 

    void prepare_differentials_kernel(mutex_type::scoped_lock& l); 

    // NOTE: Operations on each axis overlap each other.
    void compute_flux_kernel(mutex_type::scoped_lock& l);

    // NOTE: Reads from U_, writes to FX_.
    void compute_x_flux_kernel();

    // NOTE: Reads from U_, writes to FY_.
    void compute_y_flux_kernel();

    // NOTE: Reads from U_, writes to FZ_.
    void compute_z_flux_kernel();

    // IMPLEMENT 
    void adjust_flux_kernel(mutex_type::scoped_lock& l);

    void sum_differentials_kernel(mutex_type::scoped_lock& l);

  public:
    ///////////////////////////////////////////////////////////////////////////
    // IMPLEMENT
    octree_client clone_and_refine();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                clone_and_refine,
                                clone_and_refine_action);  

    ///////////////////////////////////////////////////////////////////////////
    void output();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                output,
                                output_action);  
};

}

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION_DECLARATION(                                        \
        octopus::octree_server::BOOST_PP_CAT(name, _action),                \
        BOOST_PP_CAT(octopus_octree_server_, BOOST_PP_CAT(name, _action)))  \
    /**/

// FIXME: Make sure this is in order.
OCTOPUS_REGISTER_ACTION(create_child);
OCTOPUS_REGISTER_ACTION(set_sibling);
OCTOPUS_REGISTER_ACTION(tie_sibling);
OCTOPUS_REGISTER_ACTION(set_child_sibling);
OCTOPUS_REGISTER_ACTION(tie_child_sibling);
OCTOPUS_REGISTER_ACTION(get_siblings);
OCTOPUS_REGISTER_ACTION(get_offset);
OCTOPUS_REGISTER_ACTION(get_step); // Only used by output()
OCTOPUS_REGISTER_ACTION(receive_ghost_zones);
OCTOPUS_REGISTER_ACTION(send_ghost_zone);
OCTOPUS_REGISTER_ACTION(send_mapped_ghost_zone);
OCTOPUS_REGISTER_ACTION(apply);
OCTOPUS_REGISTER_ACTION(step);
OCTOPUS_REGISTER_ACTION(step_to_time);
OCTOPUS_REGISTER_ACTION(clone_and_refine);
OCTOPUS_REGISTER_ACTION(output);

#undef OCTOPUS_REGISTER_ACTION

#endif // OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

