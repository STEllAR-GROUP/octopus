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

#include <octopus/octree/octree_client.hpp>
#include <octopus/array1d.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT octree_server
  : hpx::components::managed_component_base<octree_server>
{
  private:
    boost::array<octree_client, 8> children_;
    boost::array<octree_client, 6> siblings_;
    boost::uint64_t level_;
    array1d<boost::uint64_t, 3> location_;

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & children_;
        ar & siblings_;
        ar & level_;
        ar & location_;
    }

    void destroy_child(child_index idx);

  public:
    octree_server() : children_(), siblings_(), level_(), location_() {} 

    ~octree_server();

    ///////////////////////////////////////////////////////////////////////////
    void create_child(child_index idx);

    HPX_DEFINE_COMPONENT_ACTION(octree_server,
                                create_child,
                                create_child_action);
};

}

HPX_REGISTER_ACTION_DECLARATION(
    octopus::octree_server::create_child_action,
    octopus_octree_server_create_child_action);

#endif // OCTOPUS_58B04A8F_72F9_4B01_A8B3_941867802BA0

