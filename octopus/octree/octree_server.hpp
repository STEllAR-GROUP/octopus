////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
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
    boost::array<octree_client, 8> children_;
    boost::array<octree_client, 6> siblings_; 
    boost::uint64_t level_;
    array1d<boost::int64_t, 3> location_; ///< NOTE: Not sure if this needs to
                                          ///  be signed.

    ///////////////////////////////////////////////////////////////////////////
    // From Grid/GridNode
    double dx_;   ///< The spatial size of the node (w/o ghost zones).
                  ///  NOTE: Confirmation needed from Dominic. 
    double time_; ///< The current (physics?) time.
                  ///  NOTE: Confirmation needed from Dominic.

    array1d<boost::int64_t, 3> offset_; ///< NOTE: Not sure if this needs to be
                                        ///  signed. 

    array1d<double, 3> origin_; /// The origin of the cartesian grid. NOTE:
                                /// Confirmation needed from Dominic.

    vector3d<std::vector<double> > U_; ///< 3d array of state vectors, includes
                                       ///  ghost zones. Size of the state
                                       ///  vectors comes from the science
                                       ///  table.

    ///////////////////////////////////////////////////////////////////////////
    // TODO: Migration.
#if 0
    friend class boost::serialization::access;

    // FIXME: Should not be able to move if not initialized.
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & children_;
        ar & siblings_;
        ar & level_;
        ar & location_;
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
        if ((++siblings_set_ == 6) && received_state_)
            initialized_.set(); 
    }  

    // Precondition: mtx_ must be locked.
    void state_received_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "mutex is not locked");
        OCTOPUS_ASSERT_MSG(siblings_set_ < 6, "double initialization");
        received_state_ = true;
        if (siblings_set_ == 6)
            initialized_.set(); 
    }  

    child_index get_child_index() const
    {
        mutex_type::scoped_lock l(mtx_);
        return get_child_index_locked(l); 
    }

  public:
    octree_server()
    {
        OCTOPUS_ASSERT_MSG(false, "octree_server can't be default constructed");
    } 

    /// \brief Construct a child node.
    octree_server(
        octree_init_data const& init
      , bool root
        )
      : initialized_()
      , mtx_()
      , siblings_set_(root ? 6 : 0)
      , children_()
      , siblings_()
      , level_(init.level)
      , location_(init.location)
      , dx_(init.dx)
      , time_(init.time)
      , offset_(init.offset)
      , origin_(init.origin)
    {}

    /// \brief Construct a child node.
    octree_server(
        BOOST_RV_REF(octree_init_data) init
      , bool root
        )
      : initialized_()
      , mtx_()
      , siblings_set_(root ? 6 : 0)
      , children_()
      , siblings_()
      , level_(init.level)
      , location_(init.location)
      , dx_(init.dx)
      , time_(init.time)
      , offset_(init.offset)
      , origin_(init.origin)
    {}

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Create the \a kid child for this node.
    /// 
    /// Remote Operations:   Possibly.
    /// Concurrency Control: Waits on initialization_, locks mtx_
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
    /// Concurrency Control: Locks mtx_
    /// Synchrony Gurantee:  Fire-and-Forget 
    void set_sibling(
        boost::uint8_t f
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
        boost::uint8_t target_f
      , octree_client target_sib
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
      , boost::uint8_t f
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
      , boost::uint8_t target_f
      , octree_client target_sib
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                tie_child_sibling,
                                tie_child_sibling_action);
};

}

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::create_child_action,
    octopus_octree_server_create_child_action);

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::set_sibling_action,
    octopus_octree_server_set_sibling_action);

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::tie_sibling_action,
    octopus_octree_server_tie_sibling_action);

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::set_child_sibling_action,
    octopus_octree_server_set_child_sibling_action);

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::tie_child_sibling_action,
    octopus_octree_server_tie_child_sibling_action);

#endif // OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

