////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE)
#define OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE

#include <hpx/runtime/components/server/simple_component_base.hpp>

#include <octopus/assert.hpp>
#include <octopus/engine/configuration.hpp>
#include <octopus/engine/node_distributor_base.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT engine_server
  : hpx::components::simple_component_base<engine_server>
{
  private:
    configuration config_;

    node_distributor_base* distributor_;

  public:
    engine_server();

    bool prepare_modules(); 

    HPX_DEFINE_COMPONENT_ACTION(engine_server,
                                prepare_modules,
                                prepare_modules_action);
};

extern OCTOPUS_EXPORT engine_server* engine_ptr;

}

HPX_REGISTER_ACTION_DECLARATION(
    octopus::engine_server::prepare_modules_action,
    octopus_engine_server_prepare_modules_action);

#endif // OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE

