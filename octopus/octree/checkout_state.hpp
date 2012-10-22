////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_0D6322F2_F315_4906_BEEB_666E993F87DE)
#define OCTOPUS_0D6322F2_F315_4906_BEEB_666E993F87DE

#include <octopus/octree/octree_server.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT checkout_for_init_tag {};

extern OCTOPUS_EXPORT checkout_for_init_tag checkout_for_init;

/// This object locks its subject (an octree_server) on construction and unlocks
/// it on destruction. It also allocates the state on construction and marks the
/// node as initialized on destruction if the checkout_for_init tag is used.
/// Additionally, it provides an interface to the individual zones in the 
/// octree_server.
struct OCTOPUS_EXPORT checkout_state
{
  private:
    octree_server::mutex_type::scoped_lock lock_;
    octopus::octree_server& e_;
    bool mark_state_received_; 

  public:
    checkout_state(
        octopus::octree_server& e
        )
      : lock_(e.mtx_)
      , e_(e)
      , mark_state_received_(false)
    {}

    checkout_state(
        octopus::octree_server& e
      , checkout_for_init_tag
        )
      : lock_(e.mtx_)
      , e_(e)
      , mark_state_received_(true)
    {
        OCTOPUS_ASSERT(!e_.state_received_);
        OCTOPUS_ASSERT(e_.U_.size() == 0);
    }

    ~checkout_state()
    {
        if (mark_state_received_)
            e_.state_received_locked(lock_);
    }

    std::vector<double>& operator()(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        )
    {
        return e_.U_(i, j, k);
    }

    std::vector<double> const& operator()(
        boost::uint64_t i
      , boost::uint64_t j
      , boost::uint64_t k
        ) const
    {
        return e_.U_(i, j, k);
    }

    boost::uint64_t x_length() const
    {
        return e_.U_.x_length();
    }

    boost::uint64_t y_length() const
    {
        return e_.U_.y_length();
    }

    boost::uint64_t z_length() const
    {
        return e_.U_.z_length();
    }
};

}

#endif // OCTOPUS_0D6322F2_F315_4906_BEEB_666E993F87DE

