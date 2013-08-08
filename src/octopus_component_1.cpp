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

#undef OCTOPUS_REGISTER_ACTION

///////////////////////////////////////////////////////////////////////////////
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(
    hpx::components::simple_component<octopus::engine_server>,
    octopus_engine_server);

