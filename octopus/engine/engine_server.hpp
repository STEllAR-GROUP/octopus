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
#include <octopus/engine/science_table.hpp>
#include <octopus/array1d.hpp>
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
    config_data const config_;
    science_table const science_;

    // TODO: Replace with modular distribution.
    boost::atomic<boost::uint64_t> round_robin_;
    std::vector<hpx::id_type> const localities_;

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

    config_data const& config() const
    {
        return config_;
    }

    science_table const& science() const
    {
        return science_;
    }

    hpx::future<hpx::id_type, hpx::naming::gid_type> create_octree_async(
        octree_init_data const& init
        );

    hpx::future<hpx::id_type, hpx::naming::gid_type> create_octree_async(
        BOOST_RV_REF(octree_init_data) init
        );

    hpx::id_type create_octree(
        octree_init_data const& init
        )
    {
        return create_octree_async(init).get();
    } 

    hpx::id_type create_octree(
        BOOST_RV_REF(octree_init_data) init
        )
    {
        return create_octree_async(init).get();
    } 
};

}

#endif // OCTOPUS_2FA40BB3_E2EB_4359_842E_EF713B0FE5AE

