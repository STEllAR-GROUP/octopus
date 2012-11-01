////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_6ABA967A_917A_49FD_A3F9_3F9B69826F9F)
#define OCTOPUS_6ABA967A_917A_49FD_A3F9_3F9B69826F9F

#include <octopus/assert.hpp>

#include <hpx/runtime/components/server/simple_component_base.hpp>

#include <boost/serialization/map.hpp>
#include <boost/process.hpp>

namespace octopus
{

// Not thread-safe.
struct OCTOPUS_EXPORT visit_simulation_server
  : hpx::components::simple_component_base<visit_simulation_server>
{
  private:
    boost::process::child sim_;
    bool started_;
 
  public:
    visit_simulation_server() sim_(), started_(false) {}

    // NOTE: exec is set up by the client, as is args if none are specified.
    void start(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exec
      , std::vector<std::string> const& args
        // environment is a std::map<std::string, std::string>
      , boost::process::environment env // by value for swap
        );

    HPX_DEFINE_COMPONENT_ACTION(visit_simulation_server,
                                start,
                                start_action);

    // FIXME: Would be sweet if this could return a string with any errors from
    // the interpreter.
    void evaluate(std::string const& source); 

    HPX_DEFINE_COMPONENT_ACTION(visit_simulation_server,
                                evaluate,
                                evaluate_action);

    void terminate(); // aka crash messily

    HPX_DEFINE_COMPONENT_ACTION(visit_simulation_server,
                                terminate,
                                terminate_action);
};

}

#define OCTOPUS_REGISTER_ACTION(name)                                       \
    HPX_REGISTER_ACTION_DECLARATION(                                        \
        octopus::visit_simulation_server::BOOST_PP_CAT(name, _action),      \
        BOOST_PP_CAT(octopus_visit_simulation_server_,                      \
            BOOST_PP_CAT(name, _action)))                                   \
    /**/

OCTOPUS_REGISTER_ACTION(start);
OCTOPUS_REGISTER_ACTION(evaluate);
OCTOPUS_REGISTER_ACTION(terminate);

#undef OCTOPUS_REGISTER_ACTION

#endif // OCTOPUS_6ABA967A_917A_49FD_A3F9_3F9B69826F9F

