////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE)
#define OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE

#include <hpx/include/future.hpp>
#include <hpx/runtime/components/server/simple_component_base.hpp>

#include <octopus/octree/octree_init_data.hpp>
#include <octopus/engine/runtime_config.hpp>
#include <octopus/science/science_table.hpp>
#include <octopus/assert.hpp>

#include <boost/atomic.hpp>

#include <iostream>

// TODO: Add I/O abstraction and services.

namespace octopus
{

struct OCTOPUS_EXPORT engine_server;

extern OCTOPUS_EXPORT engine_server* engine_ptr;

struct OCTOPUS_EXPORT engine_server
  : hpx::components::simple_component_base<engine_server>
{
  private:
    config_data config_;
    science_table science_;

    // TODO: Replace with modular distribution.
    boost::atomic<boost::uint64_t> round_robin_;
    std::vector<hpx::id_type> localities_;

  public:
    engine_server() : config_(), science_(), round_robin_(0), localities_()
    {
        OCTOPUS_ASSERT_MSG(false, "engine_server can't be default constructed");
    }

    engine_server(
        config_data const& config
      , science_table const& science 
      , std::vector<hpx::id_type> const& localities
        )
      : config_(config)
      , science_(science) 
      , round_robin_(0)
      , localities_(localities)
    {
        OCTOPUS_ASSERT_MSG(engine_ptr == 0, "engine_ptr has already been set");
        engine_ptr = this;
    }

    config_data& config() 
    {
        return config_;
    }

    science_table& science() 
    {
        return science_;
    }

    std::vector<hpx::id_type> const& localities() const
    {
        return localities_;
    }

    hpx::future<hpx::id_type> create_octree_async(
        octree_init_data const& init
      , boost::shared_ptr<vector4d<double> > const& parent_U
        );

    std::vector<hpx::future<void> > call_everywhere(
        hpx::util::function<void()> const& f
        ) const;
};

}

#endif // OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE

