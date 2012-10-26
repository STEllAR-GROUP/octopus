////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_CAB36801_B41D_4CED_A034_0BA437666DB3)
#define OCTOPUS_CAB36801_B41D_4CED_A034_0BA437666DB3

#include <hpx/lcos/future_wait.hpp>
#include <hpx/async.hpp>

#include <octopus/engine/engine_interface.hpp>
#include <octopus/octree/octree_server.hpp>

namespace octopus
{

template <typename T>
inline T octree_server::reduce(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
    ) 
{
    // Make sure that we are initialized.
    initialized_.wait();

    std::vector<hpx::future<T> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    // Start recursively executing the ourself on our children.
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
            recursion_is_parallelism.push_back
                (children_[i].template reduce_async<T>(f, reducer)); 

    T result = f(*this);

    // Reduce the results from our children.
    hpx::wait(recursion_is_parallelism,
        boost::bind(&octree_server::template add_reduce<T>
                  , this, boost::ref(result), reducer, _1, _2)); 

    return result;
}

template <typename T>
inline T octree_server::reduce_zonal(
    hpx::util::function<T(std::vector<double>&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial
    ) 
{
    // Make sure that we are initialized.
    initialized_.wait();

    std::vector<hpx::future<T> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    // Start recursively executing the ourself on our children.
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
            recursion_is_parallelism.push_back
                (children_[i].template reduce_zonal_async<T>
                    (f, reducer, initial)); 

    boost::uint64_t bw = science().ghost_zone_width;
    boost::uint64_t gnx = config().grid_node_length;

    T result = initial;

    for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
        for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
            for (boost::uint64_t k = bw; k < (gnx - bw); ++k)
                result = reducer(result, f(U_(i, j, k))); 

    // Reduce the results from our children.
    hpx::wait(recursion_is_parallelism,
        boost::bind(&octree_server::template add_reduce<T>
                  , this, boost::ref(result), reducer, _1, _2)); 

    return result;
}

template <typename T>
inline T octree_client::reduce(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
    ) 
{
    return reduce_async<T>(f, reducer).get();
}

template <typename T>
inline hpx::future<T> octree_client::reduce_async(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
    ) 
{
    ensure_real();
    typedef octopus::octree_server::reduce_action<T> action_type;
    return hpx::async<action_type>(gid_, f, reducer); 
}

template <typename T>
inline T octree_client::reduce_zonal(
    hpx::util::function<T(std::vector<double>&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial 
    ) 
{
    return reduce_zonal_async<T>(f, reducer, initial).get();
}

template <typename T>
inline hpx::future<T> octree_client::reduce_zonal_async(
    hpx::util::function<T(std::vector<double>&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial
    ) 
{
    ensure_real();
    typedef octopus::octree_server::reduce_zonal_action<T> action_type;
    return hpx::async<action_type>(gid_, f, reducer, initial); 
}

}

HPX_REGISTER_ACTION_DECLARATION_TEMPLATE(
    (template <typename T>),
    (octopus::octree_server::reduce_action<T>))

HPX_REGISTER_ACTION_DECLARATION_TEMPLATE(
    (template <typename T>),
    (octopus::octree_server::reduce_zonal_action<T>))

#endif // OCTOPUS_CAB36801_B41D_4CED_A034_0BA437666DB3

