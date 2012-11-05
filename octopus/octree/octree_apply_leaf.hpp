////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_3ECF87C0_9B92_4637_A0C1_ABADBAA8CDF9)
#define OCTOPUS_3ECF87C0_9B92_4637_A0C1_ABADBAA8CDF9

#include <hpx/lcos/future_wait.hpp>
#include <hpx/async.hpp>

#include <octopus/engine/engine_interface.hpp>
#include <octopus/octree/octree_server.hpp>

// FIXME: Use boost::result_of.

namespace octopus
{

template <typename T>
inline T octree_server::apply_leaf(
    hpx::util::function<T(octree_server&)> const& f
    ) 
{
    // Make sure that we are initialized.
    //initialized_.wait();

    return f(*this);
}

template <typename T>
inline T octree_client::apply_leaf(
    hpx::util::function<T(octree_server&)> const& f
    ) const 
{
    return apply_leaf_async<T>(f).get();
}

template <typename T>
inline hpx::future<T> octree_client::apply_leaf_async(
    hpx::util::function<T(octree_server&)> const& f
    ) const 
{
    ensure_real();
    typedef octopus::octree_server::apply_leaf_action<T> action_type;
    return hpx::async<action_type>(gid_, f); 
}

}

HPX_REGISTER_ACTION_DECLARATION_TEMPLATE(
    (template <typename T>),
    (octopus::octree_server::apply_leaf_action<T>))

#endif // OCTOPUS_3ECF87C0_9B92_4637_A0C1_ABADBAA8CDF9

