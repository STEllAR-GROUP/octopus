////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/include/plain_actions.hpp>

#include <octopus/engine/engine_interface.hpp>

namespace octopus
{

void call_here(hpx::util::function<void()> const& f) 
{
    f();  
}

}

HPX_PLAIN_ACTION(octopus::call_here, call_here_action);

namespace octopus
{

std::vector<hpx::future<void> > call_everywhere(
    hpx::util::function<void()> const& f
    ) 
{
    OCTOPUS_ASSERT_MSG(!localities().empty(),
                       "no localities supporting Octopus available");

    std::vector<hpx::future<void> > calls;
    calls.reserve(localities().size());

    for (boost::uint64_t i = 0; i < localities().size(); ++i)
        calls.emplace_back(hpx::async<call_here_action>(localities()[i], f));

    return calls;
}

}


