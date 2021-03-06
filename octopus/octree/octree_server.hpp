////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bry_centere Adelstein-Lelbach
//                        ^^^^^^^ Gotta love search and replace.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0)
#define OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

#include <hpx/traits.hpp>
#include <hpx/runtime/components/server/managed_component_base.hpp>
#include <hpx/lcos/local/mutex.hpp>
#include <hpx/lcos/local/channel.hpp>

#include <octopus/array.hpp>
#include <octopus/octree/octree_init_data.hpp>
#include <octopus/octree/octree_client.hpp>
#include <octopus/atomic_bitset.hpp>

#include <bitset>

// TODO: apply_criteria, apply_zonal, apply_zonal_leaf, reduce_leaf,
// reduce_zonal_leaf
// TODO: Get rid of unnecessary _kernel and _locked suffixes.
 
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

struct OCTOPUS_EXPORT oid_type
{
  private:
    boost::uint64_t level_;
    array<boost::uint64_t, 3> location_;
    hpx::id_type gid_;

    friend std::ostream& operator<<(std::ostream& os, oid_type const& id);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & level_;
        ar & location_;
        ar & gid_; 
    }

  public:
    oid_type() : level_(0), gid_(hpx::invalid_id)
    {
        location_[0] = 0;
        location_[1] = 0;
        location_[2] = 0;
    }

    oid_type(oid_type const& other)
      : level_(other.level_), location_(other.location_), gid_(other.gid_)
    {}

    oid_type(octree_server const& e);

    oid_type& operator=(oid_type const& other)
    {
        level_ = other.level_;
        location_ = other.location_;
        gid_ = other.gid_;
        return *this;
    }

    oid_type& operator=(octree_server const& e);
};

std::ostream& operator<<(std::ostream& os, oid_type const& id);

struct OCTOPUS_EXPORT state_interpolation_data
{
    octree_client subject;
    face direction;
    array<boost::int64_t, 3> offset;

    state_interpolation_data() : subject(), direction()
    {
        offset[0] = 0;
        offset[1] = 1;
        offset[2] = 2;
    }

    state_interpolation_data(state_interpolation_data const& other)
      : subject(other.subject), direction(other.direction), offset(other.offset)
    {}

    state_interpolation_data(
        octree_client const& s
      , face d
      , array<boost::int64_t, 3> o
        )
      : subject(s), direction(d), offset(o)
    {}

    state_interpolation_data(
        octree_client const& s
      , face d
        )
      : subject(s), direction(d), offset()
    {
        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;
    }

    bool operator<(state_interpolation_data const& rhs) const
    {
        using hpx::naming::detail::strip_credit_from_gid;
        hpx::naming::gid_type lhs_gid = subject.gid_.get_gid();
        hpx::naming::gid_type rhs_gid = rhs.subject.gid_.get_gid();
        return std::make_pair(strip_credit_from_gid(lhs_gid)
                            , direction)
             < std::make_pair(strip_credit_from_gid(rhs_gid)
                            , rhs.direction);
    }

    bool operator==(state_interpolation_data const& rhs) const
    {
        using hpx::naming::detail::strip_credit_from_gid;
        hpx::naming::gid_type lhs_gid = subject.gid_.get_gid();
        hpx::naming::gid_type rhs_gid = rhs.subject.gid_.get_gid();
        return std::make_pair(strip_credit_from_gid(lhs_gid)
                            , direction)
            == std::make_pair(strip_credit_from_gid(rhs_gid)
                            , rhs.direction);
    }
};

struct OCTOPUS_EXPORT flux_interpolation_data
{
    octree_client subject;
    face direction;
    child_index idx;

    flux_interpolation_data()
      : subject(), direction(), idx()
    {}

    flux_interpolation_data(flux_interpolation_data const& other)
      : subject(other.subject), direction(other.direction), idx(other.idx) 
    {}

    flux_interpolation_data(face d, child_index i)
      : subject(), direction(d), idx(i)
    {}

    flux_interpolation_data(octree_client const& s, face d, child_index i)
      : subject(s), direction(d), idx(i) 
    {}

    bool operator<(flux_interpolation_data const& rhs) const
    {
        return std::make_pair(idx, direction)
             < std::make_pair(rhs.idx, rhs.direction); 
    }

    bool operator==(flux_interpolation_data const& rhs) const
    {
        return std::make_pair(idx, direction)
            == std::make_pair(rhs.idx, rhs.direction); 
    }
};

struct OCTOPUS_EXPORT octree_server
  : hpx::components::managed_component_base<octree_server>
{
  private:
    friend struct oid_type;

    typedef hpx::components::managed_component_base<octree_server> base_type;

    typedef hpx::components::managed_component<octree_server>*
        back_pointer_type;

    typedef hpx::lcos::local::mutex mutex_type;

    ///< This event is triggered when we are first initialized. 
//    hpx::lcos::local::event initialized_;

    mutable mutex_type mtx_; 
    hpx::id_type this_;
//    boost::uint8_t siblings_set_;
//    bool state_received_;

    // Circular doubly-linked list; size == temporal prediction gap.
    octree_client future_self_;
    octree_client past_self_;
 
//    atomic_bitset<8> marked_for_refinement_;
    std::bitset<8> marked_for_refinement_;
 
    typedef array<
        hpx::lcos::local::channel<vector4d<double> >, 6
    > sibling_state_dependencies;
  
    typedef array<
        hpx::lcos::local::channel<vector4d<double> >, 8
    > children_state_dependencies;

    typedef array<
        hpx::lcos::local::channel<vector4d<double> >, 36
    > children_flux_dependencies;

    typedef array<
        hpx::lcos::local::channel<void>, 6
    > sibling_sync_dependencies;

    // IMPLEMENT: This should totally be in the science table, along with like
    // 3k other lines of stuff in octree_server.
    /// Bryce's math for the # of communications per step (for TVD RK):
    ///
    ///     * 1 ghost zone communication at the end of each step.
    ///     * 1 ghost zone communication, 1 child -> parent state injection and 1
    ///       child -> parent flux injection during each sub step.
    ///
    /// RK1, 2 GZ comms + 1 c->p flux + 1 c->p state = 4 comms
    /// RK2, 3 GZ comms + 2 c->p flux + 2 c->p state = 7 comms 
    /// RK3, 4 GZ comms + 3 c->p flux + 3 c->p state = 10 comms 

    // Queue for incoming ghost zones.
    // NOTE: Elements of this queue should be cleared but not removed until the
    // end of each timestep. This is necessary to ensure that indices into this
    // vector remain the same throughout the entire step. 
    std::vector<sibling_state_dependencies> ghost_zone_deps_;

    // Queue for incoming state from our children.
    // NOTE: Elements of this queue should be cleared but not removed until the
    // end of each timestep. This is necessary to ensure that indices into this
    // vector remain the same throughout the entire step. 
    // NOTE: Some of these may be empty; if they are, these indicates that
    // we don't have that child.
    std::vector<children_state_dependencies> children_state_deps_;

    // Queue for incoming flux from our children.
    // NOTE: Elements of this queue should be cleared but not removed until the
    // end of each timestep. This is necessary to ensure that indices into this
    // vector remain the same throughout the entire step. 
    // NOTE: Some of these may be empty; if they are, these indicates that
    // we don't have that child.
    std::vector<children_flux_dependencies> children_flux_deps_;

    std::vector<sibling_sync_dependencies> refinement_deps_;

    ///////////////////////////////////////////////////////////////////////////
    // From OctNode
    octree_client parent_; 
    array<octree_client, 8> children_;
    array<octree_client, 6> siblings_; // FIXME: Misleading, should be
                                       // neighbors.
    std::set<state_interpolation_data> nephews_;
    std::set<flux_interpolation_data> exterior_nephews_;
    boost::uint64_t level_;
    array<boost::uint64_t, 3> location_; 

    ///////////////////////////////////////////////////////////////////////////
    // From Grid/GridNode (mostly)
    // NOTE: Rename dx_;
    double dx_;   ///< The spatial size of the node (w/o ghost zones). h in the
                  ///  original code. NOTE: Confirmation needed from Dominic.  

    double dx0_;

    // NOTE: This is only thing that keeps timesteps from overwriting each
    // other. Computing the CFL is an implicit global barrier, so I figured I
    // might as well utilize this point of synchronization to reduce memory
    // usage. P.S., only used by the root of each timestep currently. 
    hpx::lcos::local::channel<double> dt_;
    
    double time_; ///< The current (physics?) time.
                  ///  NOTE: Confirmation needed from Dominic.

    array<boost::int64_t, 3> offset_; ///< NOTE: Not sure if this needs
                                      ///  be signed. 

    array<double, 3> origin_; ///< The origin of the cartesian grid.
                                     ///  NOTE: Confirmation needed from
                                     ///  from Dominic.

    // TODO: Rename step_.
    boost::uint64_t step_;

    // REVIEW: Consider compile-time maximum sizes for the state vector, to
    // optimize allocations.
    // REVIEW: Are the "corners" of our cube are needed? I think there may be
    // eight chunks of unused space that we don't ever use.
    // 3d array of state vectors, includes ghost zones. Size of the state
    // vectors comes from the science table.
    boost::shared_ptr<vector4d<double> > U_;

    // Data from previous timestep.
    boost::shared_ptr<vector4d<double> > U0_; 

    // Scratch space for computations.
    vector4d<double> FX_; ///< Flux (X-axis).
    vector4d<double> FY_; ///< Flux (Y-axis).
    vector4d<double> FZ_; ///< Flux (Z-axis).

    boost::shared_ptr<state> FO_; ///< Flow off (stuff that
                                  ///  leaves the problem space).

    boost::shared_ptr<state> FO0_;

    // Scratch space for computations.
    vector4d<double> D_; ///< Flux differential.

    // Scratch space for computations.
    state DFO_; ///< Flow off differential. 

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
    child_index get_child_index_locked(/*mutex_type::scoped_lock& l*/) const
    {
        //OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(0 != level_, "root octree_server has no parent");
        child_index idx(location_[0] % 2, location_[1] % 2, location_[2] % 2);
        return idx; 
    }

    boost::uint8_t get_flux_index(
        octopus::axis a
      , boost::uint8_t l
      , boost::uint8_t cj
      , boost::uint8_t ck
        )
    {
        return a + l * 3 + cj * 3 * 3 + ck * 3 * 3 * 2; 
    }

    // Preconditions: mtx_ must be locked, siblings_set_ must be less than 6.
/*
    void sibling_set_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(siblings_set_ < 6, "double initialization");

        std::cout << ( boost::format("%1%: sibling_set_locked, siblings_set == %2%\n")
                     % get_gid() % boost::uint16_t(siblings_set_ + 1)); 

        if ((++siblings_set_ == 6) && state_received_)
        {
            for (std::size_t i = 0; i < 6; ++i)
                OCTOPUS_ASSERT(siblings_[i] != hpx::invalid_id);
            initialized_.set(); 
        }
    }  
*/

    // Preconditions: mtx_ must be locked, state_received_ is false. 
/*
    void state_received_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(state_received_ == false, "double initialization");

        state_received_ = true;

        std::cout << ( boost::format("%1%: sibling_set_locked, siblings_set == %2%\n")
                     % get_gid() % boost::uint16_t(siblings_set_ + 1)); 

        if (siblings_set_ == 6)
        {
            for (std::size_t i = 0; i < 6; ++i)
                OCTOPUS_ASSERT(siblings_[i] != hpx::invalid_id);
            initialized_.set();
        }
    }  
*/

    child_index get_child_index() const
    {
        //mutex_type::scoped_lock l(mtx_);
        return get_child_index_locked(/*l*/); 
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set the state of this node based on the state of its parent. 
    /// 
    /// Remote Operations:   No.
    /// Concurrency Control: Locks mtx_.
    /// Synchrony Gurantee:  Synchronous. 
    void parent_to_child_injection(
        vector4d<double> const& pU
        );

    void initialize_queues();

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Get a reference to this node that is safe to pass to our
    ///        children. Said reference must be uncounted to prevent reference
    ///        cy_centerles.
    /// 
    /// Remote Operations:   No.
    /// Concurrency Control: No.
    /// Synchrony Gurantee:  Synchronous. 
    hpx::id_type reference_from_this() const
    {
        // We shouldn't need to lock here, I believe.
//        hpx::id_type gid = get_gid();
        OCTOPUS_ASSERT(hpx::invalid_id != this_);
        OCTOPUS_ASSERT_MSG(
            this_.get_management_type() == hpx::id_type::unmanaged,
            "get_gid() should return an unmanaged GID");
        return this_;
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Get a client referencing this node that is safe to pass to our
    ///        children. The reference the child holds must be uncounted to
    ///        prevent reference cy_centerles.
    /// 
    /// Remote Operations:   No.
    /// Concurrency Control: No.
    /// Synchrony Gurantee:  Synchronous. 
    octree_client client_from_this() const
    {
        // We shouldn't need to lock here, I believe.
//        hpx::id_type gid = get_gid();
        OCTOPUS_ASSERT(hpx::invalid_id != this_);
        OCTOPUS_ASSERT_MSG(
            this_.get_management_type() == hpx::id_type::unmanaged,
            "get_gid() should return an unmanaged GID");
        return octree_client(this_);
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

    /// \brief Construct a child node.
    octree_server(
        back_pointer_type back_ptr
      , octree_init_data const& init
      , boost::shared_ptr<vector4d<double> > const& parent_U
        );

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

    double get_dx() const
    {
        return dx_;
    }

    double get_dt() const
    {
        return dt_.get();
    }

    void post_dt(double dt) 
    {
        OCTOPUS_ASSERT(0 == level_);
        dt_.reset();
        dt_.post(dt);
    }

    state& operator()(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        )
    {
        return (*U_)(i, j, k);
    } 

    state const& operator()(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
        return (*U_)(i, j, k);
    } 

/*
    // NOTE: Use with caution.
    mutex_type& get_mutex() const
    {
        return mtx_;
    }
*/

    array<double, 3> center_coords(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
        array<double, 3> coords;
        coords[0] = x_center(i);
        coords[1] = y_center(j);
        coords[2] = z_center(k);
        return coords;
    }

    array<double, 3> x_face_coords(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
    	array<double, 3> coords;
        coords[0] = x_face(i);
        coords[1] = y_center(j);
        coords[2] = z_center(k);
    	return coords;
    }

    array<double, 3> y_face_coords(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
    	array<double, 3> coords;
        coords[0] = x_center(i);
        coords[1] = y_face(j);
        coords[2] = z_center(k);
    	return coords;
    }
    
    array<double, 3> z_face_coords(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
    	array<double, 3> coords;
        coords[0] = x_center(i);
        coords[1] = y_center(j);
        coords[2] = z_face(k);
    	return coords;
    }

    double x_center(boost::uint64_t i) const
    {
        return x_face(i) + 0.5 * dx_;
    }

    double y_center(boost::uint64_t j) const
    {
        return y_face(j) + 0.5 * dx_;
    }

    double z_center(boost::uint64_t k) const
    {
        return z_face(k) + 0.5 * dx_;
    }

    double x_face(boost::uint64_t i) const;

    double y_face(boost::uint64_t i) const;

    double z_face(boost::uint64_t i) const;

    std::ostream& debug() const
    {
        return std::cout << get_oid() << ": ";
    }

    ///////////////////////////////////////////////////////////////////////////
    void prepare_compute_queues();

    void set_time(
        double time
      , boost::uint64_t step
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                set_time,
                                set_time_action);

    void set_buffer_links(
        hpx::id_type const& future_self
      , hpx::id_type const& past_self
        )
    {
        future_self_ = octree_client(future_self);
        past_self_ = octree_client(past_self);
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                set_buffer_links,
                                set_buffer_links_action);

    void clear_refinement_marks();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                clear_refinement_marks,
                                clear_refinement_marks_action);

  private:
    void clear_refinement_marks_kernel();

  public:
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

    void require_child(
        child_index kid
        )
    {
        //debug() << "require_child(" << kid << ")\n";
        // REVIEW: Lock here.
        mutex_type::scoped_lock l(mtx_);

        if (hpx::invalid_id == children_[kid])
        {
            marked_for_refinement_.set(kid, true);
            propagate_locked(kid, l);
        }
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                require_child,
                                require_child_action);

    // FIXME: Use remotable futures.
    void require_sibling_child(
        child_index kid
      , face f
        )
    {
        if (siblings_[f].kind() == amr_boundary)
            siblings_[f].require_child(kid); 
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                require_sibling_child,
                                require_sibling_child_action);

    // FIXME: Use remotable futures.
    void require_corner_child(
        child_index kid
      , face f0
      , face f1
        )
    {
        siblings_[f0].require_sibling_child(kid, f1); 
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                require_corner_child,
                                require_corner_child_action);
    void remove_nephew(
        octree_client const& nephew
      , face f
      , child_index idx
        )
    {
        mutex_type::scoped_lock l(mtx_);
        state_interpolation_data sid(nephew.gid_, f);

        //bool erased = nephews_.erase(sid) != 0;
        //OCTOPUS_ASSERT(erased); 

        nephews_.erase(sid);
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                remove_nephew,
                                remove_nephew_action);

  public:
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
    oid_type get_oid() const
    {
        return oid_type(*this); 
    }

    HPX_DEFINE_COMPONENT_CONST_ACTION(octree_server,
                                       get_oid,
                                       get_oid_action);

    ///////////////////////////////////////////////////////////////////////////
    array<octree_client, 6> get_siblings() 
    {
        // Make sure that we are initialized.
        //initialized_.wait();

        // Is the lock needed?
        //mutex_type::scoped_lock l(mtx_);
        return siblings_;
    }

    HPX_DEFINE_COMPONENT_CONST_ACTION(octree_server,
                                get_siblings,
                                get_siblings_action);

    ///////////////////////////////////////////////////////////////////////////
    // FIXME: Remove the need for this.
    array<boost::int64_t, 3> get_offset() const
    {
        return offset_;
    }

    HPX_DEFINE_COMPONENT_CONST_ACTION(octree_server,
                                             get_offset,
                                             get_offset_action);

    ///////////////////////////////////////////////////////////////////////////
    // Purely for debugging. 
    array<boost::uint64_t, 3> get_location() const
    {
        return location_;
    }

    HPX_DEFINE_COMPONENT_CONST_ACTION(octree_server,
                                      get_location,
                                      get_location_action);

    ///////////////////////////////////////////////////////////////////////////
    // Ghost zone communication
    // NOTE: This was contained in enforce_boundaries in the original code.

    // REVIEW: I think step 2.) can come before step 1.).
    /// 0.) Push ghost zone data to our siblings and determine which ghost zones
    ///     we will receive.
    /// 1.) Wait for our ghost zones to be delivered by our siblings.
    /// 2.) Push ghost zone data to our nephews.
    void communicate_ghost_zones(
        boost::uint64_t phase
        );

  private:
    void add_ghost_zone(
        face f
      , BOOST_RV_REF(vector4d<double>) zone
        );

    void add_ghost_zone_callback(
        face f ///< Bound parameter.
      , hpx::future<vector4d<double> > zone_f
        )
    {
        add_ghost_zone(f, boost::move(zone_f.move()));
    }

  public:
    // FIXME: Rvalue reference kung-fo must be applied here.
    /// Called by our siblings.
    void receive_ghost_zone(
        boost::uint64_t step ///< For debugging purposes.
      , boost::uint64_t phase 
      , face f ///< Relative to caller.
      , vector4d<double> const& zone
        )
    {
//        mutex_type::scoped_lock l(mtx_);

        OCTOPUS_ASSERT_MSG(step_ == step,
            "cross-timestep communication occurred, octree is ill-formed");

        OCTOPUS_ASSERT_FMT_MSG(
            phase < ghost_zone_deps_.size(),
            "phase (%1%) is greater than the ghost zone queue length (%2%)",
            phase % ghost_zone_deps_.size());

        OCTOPUS_ASSERT(f != invalid_face);

        // NOTE (wash): boost::move should be safe here, zone is a temporary,
        // even if we're local to the caller. Plus, ATM set_value requires the
        // value to be moved to it.
        ghost_zone_deps_[phase](f).post(zone);
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                receive_ghost_zone,
                                receive_ghost_zone_action);

  private:
    vector4d<double> send_ghost_zone_locked(
        face f ///< Our direction, relative to the caller.
      /*, mutex_type::scoped_lock& l*/
        );
  public:

    /// Produces ghost zone data for a sibling.
    vector4d<double> send_ghost_zone(
        face f ///< Our direction, relative to the caller.
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                send_ghost_zone,
                                send_ghost_zone_action);

    // FIXME: The octree doing the interpolation could reverse-engineer the
    // offset (probably) from the face (the child_index may also be needed).
    vector4d<double> send_interpolated_ghost_zone(
        face f ///< Our direction, relative to the caller.
//      , boost::uint64_t disparity ///< Difference in refinement level
      , array<boost::int64_t, 3> offset
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                send_interpolated_ghost_zone,
                                send_interpolated_ghost_zone_action);

    void map_ghost_zone(
        face f ///< Our direction, relative to the caller.
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                map_ghost_zone,
                                map_ghost_zone_action);

    ///////////////////////////////////////////////////////////////////////////
    // Child -> parent injection of state.
    void child_to_parent_state_injection(
        boost::uint64_t phase
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                child_to_parent_state_injection,
                                child_to_parent_state_injection_action);

  private:
    /// 0.) Blocks until the child -> parent state injection that is at \a phase
    ///     in the queue is ready.
    /// 1.) Send a child -> parent state injection up to our parent. 
    void child_to_parent_state_injection_kernel(
        boost::uint64_t phase
        );

    // FIXME: Rvalue reference kung-fo must be applied here.
    /// Callback used to wait for a particular child state. 
    void add_child_state(
        child_index idx ///< Bound parameter.
      , hpx::future<vector4d<double> > state_f
        );

  public:
    // FIXME: Rvalue reference kung-fo must be applied here.
    /// Called by our children.
    void receive_child_state(
        boost::uint64_t step ///< For debugging purposes.
      , boost::uint64_t phase 
      , child_index idx 
      , BOOST_RV_REF(vector4d<double>) s
        )
    { // {{{
        //mutex_type::scoped_lock l(mtx_);

        OCTOPUS_ASSERT_MSG(step_ == step,
            "cross-timestep communication occurred, octree is ill-formed");

        OCTOPUS_ASSERT_FMT_MSG(
            phase < children_state_deps_.size(),
            "phase (%1%) is greater than the children state queue length (%2%)",
            phase % children_state_deps_.size());

        OCTOPUS_ASSERT(boost::uint64_t(idx) < 8);

        // NOTE (wash): boost::move should be safe here, zone is a temporary,
        // even if we're local to the caller. Plus, ATM set_value requires the
        // value to be moved to it.
        children_state_deps_[phase](idx).post(s);
    } // }}}

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                receive_child_state,
                                receive_child_state_action);

  private:
    vector4d<double> send_child_state();

  public:
    ///////////////////////////////////////////////////////////////////////////
    // Child -> parent injection of flux.
    void child_to_parent_flux_injection(
        boost::uint64_t phase
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                child_to_parent_flux_injection,
                                child_to_parent_flux_injection_action);

  private:
    /// 0.) Blocks until the child -> parent flux injection that is at \a phase
    ///     in the queue is ready.
    /// 1.) Send a child -> parent flux injection up to our parent. 
    void child_to_parent_flux_injection_kernel(
        boost::uint64_t phase
        );

    // FIXME: Rvalue reference kung-fo must be applied here.
    /// Callback used to wait for a particular child flux. 
    void add_child_flux(
        axis a ///< Bound parameter.
      , boost::uint64_t i0 ///< Bound parameter.
      , boost::uint8_t cj ///< Bound parameter.
      , boost::uint8_t ck ///< Bound parameter.
      , hpx::future<vector4d<double> > flux_f
        );

  public:
    // FIXME: Rvalue reference kung-fo must be applied here.
    /// Called by our children.
    void receive_child_flux(
        boost::uint64_t step ///< For debugging purposes.
      , boost::uint64_t phase 
      , boost::uint8_t idx 
      , BOOST_RV_REF(vector4d<double>) s
        )
    { // {{{
        OCTOPUS_ASSERT_MSG(step_ == step,
            "cross-timestep communication occurred, octree is ill-formed");

        OCTOPUS_ASSERT_FMT_MSG(
            phase < children_flux_deps_.size(),
            "phase (%1) is greater than the children flux queue length (%2%)",
            phase % children_flux_deps_.size());

        OCTOPUS_ASSERT(idx < 36);

        // NOTE (wash): boost::move should be safe here, zone is a temporary,
        // even if we're local to the caller. Plus, ATM set_value requires the
        // value to be moved to it.
        children_flux_deps_[phase](idx).post(s);
    } // }}}

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                receive_child_flux,
                                receive_child_flux_action);

  private:
    vector4d<double> send_child_flux(face f);

  public:
    ///////////////////////////////////////////////////////////////////////////
    void step_recurse(double dt);

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                step_recurse,
                                step_recurse_action);  

    void step()
    {
        OCTOPUS_ASSERT(0 == level_);
        step_recurse(get_dt());
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                step,
                                step_action);  

  private:
    void step_kernel(double dt);

    void sub_step_kernel(boost::uint64_t phase, double dt, double beta);

    void add_differentials_kernel(double dt, double beta); 

    void prepare_differentials_kernel(); 

    // Operations on each axis overlap each other.
    void compute_flux_kernel(boost::uint64_t phase);

    // Reads from U_, writes to FX_.
    void compute_x_flux_kernel();

    // Reads from U_, writes to FY_.
    void compute_y_flux_kernel();

    // Reads from U_, writes to FZ_.
    void compute_z_flux_kernel();

    void sum_differentials_kernel();

  public:
    ///////////////////////////////////////////////////////////////////////////
    void copy_and_regrid();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                copy_and_regrid,
                                copy_and_regrid_action);  

    void mark();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                mark,
                                mark_action);  

    void populate();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                populate,
                                populate_action);  

    void link();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                link,
                                link_action);  

    // FIXME: Remove the need for this + the second populate/link pass.
    void remark();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                remark,
                                remark_action);  

  private:
    void mark_kernel();

    void propagate_locked(child_index kid, mutex_type::scoped_lock &l);

    void populate_kernel();

    void link_kernel();

    void link_child(
        std::vector<hpx::future<void> >& links
      , child_index kid
        );

    void remark_kernel();

  public:
    void refine();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                refine,
                                refine_action);  

  private:
    void sibling_refinement_signal(boost::uint64_t phase);

  public:
    void receive_sibling_refinement_signal(boost::uint64_t phase, face f)
    {
        mutex_type::scoped_lock l(mtx_);

        OCTOPUS_ASSERT(invalid_face != f);

        refinement_deps_[phase](f).post();
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                receive_sibling_refinement_signal, 
                                receive_sibling_refinement_signal_action);

    ///////////////////////////////////////////////////////////////////////////
    void output()
    {
        client_from_this().output();
    }

    void output(double time)
    {
        client_from_this().output(time);
    }

    ///////////////////////////////////////////////////////////////////////////
    // {{{ apply - definitions are out-of-line in octree_apply.hpp.
    void apply(
        hpx::util::function<void(octree_server&)> const& f
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                apply,
                                apply_action);  
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ apply_leaf - definitions are out-of-line in octree_apply_leaf.hpp.
    template <typename T>
    T apply_leaf(
        hpx::util::function<T(octree_server&)> const& f
        );

    template <typename T>
    struct apply_leaf_action
      : hpx::actions::make_action<
            T (octree_server::*)(hpx::util::function<T(octree_server&)> const&)
          , &octree_server::template apply_leaf<T>
          , apply_leaf_action<T>
        >
    {};
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ reduce - definitions are out-of-line in octree_reduce.hpp.
  private:
    template <typename T>
    void add_reduce(
        T& result
      , hpx::util::function<T(T const&, T const&)> const& reducer
      , hpx::future<T> value
        )
    {
        T tmp = value.get();

        mutex_type::scoped_lock l(mtx_);
        result = reducer(result, boost::move(tmp)); 
    }

  public:
    // Initial should be an identity.
    template <typename T>
    T reduce(
        hpx::util::function<T(octree_server&)> const& f
      , hpx::util::function<T(T const&, T const&)> const& reducer
      , T const& initial = T()
        );

    template <typename T>
    struct reduce_action
      : hpx::actions::make_action<
            T (octree_server::*)
                ( hpx::util::function<T(octree_server&)> const&
                , hpx::util::function<T(T const&, T const&)> const&
                , T const&)
          , &octree_server::template reduce<T>
          , reduce_action<T>
        >
    {};

    // Initial should be an identity.
    template <typename T>
    T reduce_zonal(
        hpx::util::function<T(state&)> const& f
      , hpx::util::function<T(T const&, T const&)> const& reducer
      , T const& initial = T()
        );

    template <typename T>
    struct reduce_zonal_action
      : hpx::actions::make_action<
            T (octree_server::*)
                ( hpx::util::function<T(state&)> const&
                , hpx::util::function<T(T const&, T const&)> const&
                , T const&)
          , &octree_server::template reduce_zonal<T>
          , reduce_zonal_action<T>
        >
    {};
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    void slice(
        slice_function const& f
      , axis a
      , double eps = std::numeric_limits<double>::epsilon()
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                slice,
                                slice_action);

  private:
    void slice_x_kernel(slice_function const& f, double eps);

    void slice_y_kernel(slice_function const& f, double eps);

    void slice_z_kernel(slice_function const& f, double eps);

  public:
    void slice_leaf(
        slice_function const& f
      , axis a
      , double eps = std::numeric_limits<double>::epsilon()
        )
    {
        switch (a)
        {
            case axis::x_axis:
                slice_x_kernel(f, eps);
                break;
            case axis::y_axis:
                slice_y_kernel(f, eps);
                break;
            case axis::z_axis:
                slice_z_kernel(f, eps);
                break;
            default: OCTOPUS_ASSERT(false); break; 
        };
    }

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                slice_leaf,
                                slice_leaf_action);

    ///////////////////////////////////////////////////////////////////////////
    void save();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                save,
                                save_action);

    void load();

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                load,
                                load_action);
};

inline oid_type::oid_type(octree_server const& e)
  : level_(e.level_)
  , location_(e.location_)
  , gid_(e.reference_from_this())
{ }

inline oid_type& oid_type::operator=(octree_server const& e)
{
    level_ = e.level_;
    location_ = e.location_;
    gid_ = e.reference_from_this();
    return *this;
}

}

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION_DECLARATION(                                        \
        octopus::octree_server::BOOST_PP_CAT(name, _action),                \
        BOOST_PP_CAT(octopus_octree_server_, BOOST_PP_CAT(name, _action)))  \
    /**/

// FIXME: Make sure this is in order.
OCTOPUS_REGISTER_ACTION(set_time);
OCTOPUS_REGISTER_ACTION(set_buffer_links);
OCTOPUS_REGISTER_ACTION(clear_refinement_marks);

OCTOPUS_REGISTER_ACTION(create_child);
OCTOPUS_REGISTER_ACTION(require_child);
OCTOPUS_REGISTER_ACTION(require_sibling_child);
OCTOPUS_REGISTER_ACTION(require_corner_child);
OCTOPUS_REGISTER_ACTION(remove_nephew);
OCTOPUS_REGISTER_ACTION(set_sibling);
OCTOPUS_REGISTER_ACTION(tie_sibling);
OCTOPUS_REGISTER_ACTION(set_child_sibling);
OCTOPUS_REGISTER_ACTION(tie_child_sibling);

OCTOPUS_REGISTER_ACTION(get_oid);
OCTOPUS_REGISTER_ACTION(get_siblings);
OCTOPUS_REGISTER_ACTION(get_offset);
OCTOPUS_REGISTER_ACTION(get_location);

OCTOPUS_REGISTER_ACTION(receive_ghost_zone);
OCTOPUS_REGISTER_ACTION(send_ghost_zone);
OCTOPUS_REGISTER_ACTION(send_interpolated_ghost_zone);
OCTOPUS_REGISTER_ACTION(map_ghost_zone);

OCTOPUS_REGISTER_ACTION(child_to_parent_state_injection);
OCTOPUS_REGISTER_ACTION(receive_child_state);

OCTOPUS_REGISTER_ACTION(child_to_parent_flux_injection);
OCTOPUS_REGISTER_ACTION(receive_child_flux);

OCTOPUS_REGISTER_ACTION(apply);

OCTOPUS_REGISTER_ACTION(step);
OCTOPUS_REGISTER_ACTION(step_recurse);

OCTOPUS_REGISTER_ACTION(copy_and_regrid);
OCTOPUS_REGISTER_ACTION(refine);
OCTOPUS_REGISTER_ACTION(mark);
OCTOPUS_REGISTER_ACTION(populate);
OCTOPUS_REGISTER_ACTION(link);
OCTOPUS_REGISTER_ACTION(remark);
OCTOPUS_REGISTER_ACTION(receive_sibling_refinement_signal);

OCTOPUS_REGISTER_ACTION(slice);
OCTOPUS_REGISTER_ACTION(slice_leaf);

#undef OCTOPUS_REGISTER_ACTION

#endif // OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

