////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/runtime/components/stubs/runtime_support.hpp>

#include <octopus/engine/engine_server.hpp>
#include <octopus/octree/octree_server.hpp>

namespace octopus
{

engine_server* engine_ptr = 0;

hpx::future<hpx::id_type, hpx::naming::gid_type>
engine_server::create_octree_async(
    boost::uint64_t level
  , array1d<boost::uint64_t, 3> const& location
    )
{
    OCTOPUS_ASSERT_MSG(!localities_.empty(),
                       "no localities supporting Octopus available");

    using hpx::components::stubs::runtime_support;
    return runtime_support::create_component_async<octopus::octree_server>
        (localities_[round_robin_++ % localities_.size()], level, location);
}

}

