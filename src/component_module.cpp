////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/runtime/components/component_factory.hpp>

#include <hpx/util/portable_binary_iarchive.hpp>
#include <hpx/util/portable_binary_oarchive.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/export.hpp>

#include <octopus/octree/octree_server.hpp>

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(
    hpx::components::managed_component<octopus::octree_server>,
    octopus_octree_server);

HPX_REGISTER_ACTION(
    octopus::octree_server::create_child_action,
    octopus_octree_server_create_child_action);

HPX_REGISTER_ACTION(
    octopus::octree_server::set_sibling_action,
    octopus_octree_server_set_sibling_action);

HPX_REGISTER_ACTION(
    octopus::octree_server::tie_sibling_action,
    octopus_octree_server_tie_sibling_action);

HPX_REGISTER_ACTION(
    octopus::octree_server::set_child_sibling_action,
    octopus_octree_server_set_child_sibling_action);

HPX_REGISTER_ACTION(
    octopus::octree_server::tie_child_sibling_action,
    octopus_octree_server_tie_child_sibling_action);

