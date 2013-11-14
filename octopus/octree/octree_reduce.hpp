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

// Reimplement with when_all

namespace octopus
{

template <typename T>
inline T octree_server::reduce(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial
    ) 
{
    // Make sure that we are initialized.
    //initialized_.wait();

    std::vector<hpx::future<T> > keep_alive;
    keep_alive.reserve(8);

    std::vector<hpx::future<void> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    T result = initial;

    // Start recursively executing the ourself on our children.
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
        {
            keep_alive.emplace_back
                (children_[i].template reduce_async<T>(f, reducer, initial));

            // Reduce the results from our children.
            recursion_is_parallelism.emplace_back(
                keep_alive.back().then(boost::bind(
                    &octree_server::template add_reduce<T>
                  , this, boost::ref(result), reducer, _1))); 
        }

    T local_result = f(*this);

    hpx::wait(recursion_is_parallelism);

    return reducer(result, local_result);
}

// TODO: Figure out what the guranteed order actually is and document it.
template <typename T>
inline T octree_server::reduce_ordered(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial
    ) 
{
    std::vector<hpx::future<T> > keep_alive;
    keep_alive.reserve(8);

    std::vector<hpx::future<void> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    boost::array<T, 8> results;

    // Start recursively executing the ourself on our children.
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
        {
            keep_alive.emplace_back
                (children_[i].template reduce_ordered_async<T>
                    (f, reducer, initial));

            // Buffer the results from our children.
            recursion_is_parallelism.emplace_back(
                keep_alive.back().then(boost::bind(
                    &octree_server::template buffer_reduce_ordered<T>
                  , this, i, boost::ref(results), _1))); 
        }

    T local_result = f(*this);

    hpx::wait(recursion_is_parallelism);

    T result = initial;

    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
            result = reducer(result, results[i]);

    return reducer(result, local_result);
}

template <typename T>
inline T octree_server::reduce_zonal(
    hpx::util::function<T(state&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial
    ) 
{
    // Make sure that we are initialized.
    //initialized_.wait();

    std::vector<hpx::future<T> > keep_alive;
    keep_alive.reserve(8);

    std::vector<hpx::future<void> > recursion_is_parallelism;
    recursion_is_parallelism.reserve(8);

    T result = initial;

    // Start recursively executing the ourself on our children.
    for (boost::uint64_t i = 0; i < 8; ++i)
        if (hpx::invalid_id != children_[i])
        {
            keep_alive.emplace_back
                (children_[i].template reduce_zonal_async<T>
                    (f, reducer, initial)); 

            // Reduce the results from our children.
            recursion_is_parallelism.emplace_back(
                keep_alive.back().then(boost::bind(
                    &octree_server::template add_reduce<T>
                  , this, boost::ref(result), reducer, _1))); 
        }

    boost::uint64_t bw = science().ghost_zone_length;
    boost::uint64_t gnx = config().grid_node_length;

    T local_result = initial;

    for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
        for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
            for (boost::uint64_t k = bw; k < (gnx - bw); ++k)
                local_result = reducer(local_result, f((*U_)(i, j, k))); 

    hpx::wait(recursion_is_parallelism);

    return reducer(result, local_result);
}

template <typename T>
inline T octree_client::reduce(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial 
    ) const
{
    return reduce_async<T>(f, reducer, initial).get();
}

template <typename T>
inline hpx::future<T> octree_client::reduce_async(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial 
    ) const
{
    ensure_real();
    typedef octopus::octree_server::reduce_action<T> action_type;
    return hpx::async<action_type>(gid_, f, reducer, initial); 
}

template <typename T>
inline T octree_client::reduce_ordered(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial 
    ) const
{
    return reduce_ordered_async<T>(f, reducer, initial).get();
}

template <typename T>
inline hpx::future<T> octree_client::reduce_ordered_async(
    hpx::util::function<T(octree_server&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial 
    ) const
{
    ensure_real();
    typedef octopus::octree_server::reduce_ordered_action<T> action_type;
    return hpx::async<action_type>(gid_, f, reducer, initial); 
}

template <typename T>
inline T octree_client::reduce_zonal(
    hpx::util::function<T(state&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial 
    ) const 
{
    return reduce_zonal_async<T>(f, reducer, initial).get();
}

template <typename T>
inline hpx::future<T> octree_client::reduce_zonal_async(
    hpx::util::function<T(state&)> const& f
  , hpx::util::function<T(T const&, T const&)> const& reducer
  , T const& initial
    ) const
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
    (octopus::octree_server::reduce_ordered_action<T>))

HPX_REGISTER_ACTION_DECLARATION_TEMPLATE(
    (template <typename T>),
    (octopus::octree_server::reduce_zonal_action<T>))

#endif // OCTOPUS_CAB36801_B41D_4CED_A034_0BA437666DB3

