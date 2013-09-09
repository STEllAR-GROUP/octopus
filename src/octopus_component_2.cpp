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

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION(                                                    \
        octopus::octree_server::BOOST_PP_CAT(name, _action),                \
        BOOST_PP_CAT(octopus_octree_server_, BOOST_PP_CAT(name, _action)))  \
    /**/

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

OCTOPUS_REGISTER_ACTION(save);
OCTOPUS_REGISTER_ACTION(load);

#undef OCTOPUS_REGISTER_ACTION

