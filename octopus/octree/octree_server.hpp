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

    boost::array<octree_client, 8> children_;
    boost::array<octree_client, 6> siblings_;
    boost::uint64_t level_;
    array1d<boost::uint64_t, 3> location_;

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

    // Precondition: mtx_ must be locked.
    child_index get_child_index_locked(mutex_type::scoped_lock& l) const
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "node mutex is not locked");
        OCTOPUS_ASSERT_MSG(0 != level_, "root node has no parent");
        child_index idx(location_[0] % 2, location_[1] % 2, location_[2] % 2);
        return idx; 
    }

    // Precondition: mtx_ must be locked.
    void initialize_if_ready_locked(mutex_type::scoped_lock& l)
    {
        OCTOPUS_ASSERT_MSG(l.owns_lock(), "node mutex is not locked");
        OCTOPUS_ASSERT_MSG(siblings_set_ < 6, "node is already initialized");
        if (++siblings_set_ == 6)
            initialized_.set(); 
    }  

  public:
    /// \brief Construct a root node.
    octree_server()
      : initialized_()
      , mtx_()
      , siblings_set_(0)
      , children_()
      , siblings_()
      , level_()
      , location_()
    {
        initialized_.set();
    } 

    /// \brief Construct a child node.
    octree_server(
        boost::uint64_t level
      , array1d<boost::uint64_t, 3> const& location
        )
      : initialized_()
      , mtx_()
      , siblings_set_(0)
      , children_()
      , siblings_()
      , level_(level)
      , location_(location)
    {}

    child_index get_child_index() const
    {
        mutex_type::scoped_lock l(mtx_);
        return get_child_index_locked(l); 
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Create the \a kid child for this node.
    /// 
    /// Communication:       Local and possibly remote
    /// Concurrency Control: Waits on initialization_, locks mtx_
    /// Synchrony Gurantees: Fire-and-Forget 
    void create_child(
        child_index kid
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                create_child,
                                create_child_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Set \a target_sib as the \a target_f sibling of this node.
    /// 
    /// Communication:       Local and possibly remote
    /// Concurrency Control: Locks mtx_
    /// Synchrony Gurantees: Fire-and-Forget 
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
    /// Communication:       Local and possibly remote.
    /// Concurrency Control: Locks mtx_.
    /// Synchrony Gurantees: Fire-and-Forget. 
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
    /// Communication:       Local and possibly remote.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantees: Fire-and-Forget. 
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
    /// Communication:       Local and possibly remote.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantees: Fire-and-Forget.
    void tie_child_sibling(
        child_index target_kid
      , boost::uint8_t target_f
      , octree_client target_sib
        );

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                tie_child_sibling,
                                tie_child_sibling_action);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Count the number of nodes under this node. 
    /// 
    /// Communication:       Local and possibly remote.
    /// Concurrency Control: Waits on initialization_, locks mtx_.
    /// Synchrony Gurantees: Synchronous.
    boost::uint64_t get_node_count();

    HPX_DEFINE_COMPONENT_DIRECT_ACTION(octree_server,
                                       get_node_count,
                                       get_node_count_action);
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

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::get_node_count_action,
    octopus_octree_server_get_node_count_action);


#endif // OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

