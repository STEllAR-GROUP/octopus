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

#include <octopus/visit/visit_simulation_server.hpp>

HPX_REGISTER_COMPONENT_MODULE();

///////////////////////////////////////////////////////////////////////////////
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(
    hpx::components::simple_component<octopus::visit_simulation_server>,
    octopus_visit_simulation_server);

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION(                                                    \
        octopus::visit_simulation_server::BOOST_PP_CAT(name, _action),      \
        BOOST_PP_CAT(octopus_visit_simulation_server_,                      \
        BOOST_PP_CAT(name, _action)))                                       \
    /**/

OCTOPUS_REGISTER_ACTION(start);
OCTOPUS_REGISTER_ACTION(evaluate);
OCTOPUS_REGISTER_ACTION(terminate);

#undef OCTOPUS_REGISTER_ACTION

