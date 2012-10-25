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
#include <octopus/octree/octree_reduce.hpp>
#include <octopus/engine/engine_server.hpp>

HPX_REGISTER_COMPONENT_MODULE();

///////////////////////////////////////////////////////////////////////////////
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(
    hpx::components::managed_component<octopus::octree_server>,
    octopus_octree_server);

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION(                                                    \
        octopus::octree_server::BOOST_PP_CAT(name, _action),                \
        BOOST_PP_CAT(octopus_octree_server_, BOOST_PP_CAT(name, _action)))  \
    /**/

OCTOPUS_REGISTER_ACTION(create_child);
OCTOPUS_REGISTER_ACTION(set_sibling);
OCTOPUS_REGISTER_ACTION(tie_sibling);
OCTOPUS_REGISTER_ACTION(set_child_sibling);
OCTOPUS_REGISTER_ACTION(tie_child_sibling);
OCTOPUS_REGISTER_ACTION(get_siblings);
OCTOPUS_REGISTER_ACTION(get_offset);
OCTOPUS_REGISTER_ACTION(receive_ghost_zones);
OCTOPUS_REGISTER_ACTION(send_ghost_zone);
OCTOPUS_REGISTER_ACTION(send_mapped_ghost_zone);
OCTOPUS_REGISTER_ACTION(apply);
OCTOPUS_REGISTER_ACTION(apply_leaf);
OCTOPUS_REGISTER_ACTION(step);
OCTOPUS_REGISTER_ACTION(step_to_time);
OCTOPUS_REGISTER_ACTION(clone_and_refine);
OCTOPUS_REGISTER_ACTION(output);

#undef OCTOPUS_REGISTER_ACTION

///////////////////////////////////////////////////////////////////////////////
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(
    hpx::components::simple_component<octopus::engine_server>,
    octopus_engine_server);

