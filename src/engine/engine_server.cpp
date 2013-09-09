////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/runtime/components/stubs/runtime_support.hpp>
#include <hpx/include/plain_actions.hpp>

#include <octopus/engine/engine_server.hpp>
#include <octopus/octree/octree_server.hpp>
#include <octopus/trivial_serialization.hpp>
#include <octopus/science.hpp>

namespace octopus
{

engine_server* engine_ptr = 0;

engine_server::engine_server(
    config_data const& config
  , science_table const& science 
  , std::vector<hpx::id_type> const& localities
    )
  : config_(config)
  , science_(science) 
  , round_robin_(0)
  , localities_(localities)
  , checkpoint_file_()
{
    OCTOPUS_ASSERT_MSG(engine_ptr == 0, "engine_ptr has already been set");
    engine_ptr = this;

    open_checkpoint(config_.checkpoint_file, config_.load_checkpoint);
}

hpx::future<hpx::id_type> engine_server::create_octree_async(
    octree_init_data const& init
  , boost::shared_ptr<vector4d<double> > const& parent_U
    )
{
    OCTOPUS_ASSERT_MSG(!localities_.empty(),
                       "no localities supporting Octopus available");

    using hpx::components::stubs::runtime_support;

    hpx::id_type locality = science().distribute(init, localities_);

    return runtime_support::create_component_async<octopus::octree_server>
        (locality, init, parent_U);
}

void engine_server::open_checkpoint(
    std::string const& file_name
  , bool load 
    )
{
    OCTOPUS_ASSERT(!checkpoint_file_.is_open());

    std::fstream::openmode m;

    if (load)
        m = std::fstream::binary | std::fstream::in | std::fstream::in;
    else
        m = std::fstream::binary | std::fstream::out | std::fstream::trunc;

    try
    {
        std::string s = boost::str( boost::format(file_name)
                                  % hpx::get_locality_id());
        checkpoint_file_.open(s, m); 
    }
    // FIXME: Catch the specific boost.format exception.
    catch (...)
    {
        checkpoint_file_.open(file_name, m); 
    }

    OCTOPUS_ASSERT(checkpoint_file_.is_open());
}

}

