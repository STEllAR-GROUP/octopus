////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//  Copyright (c) 2007-2012 Hartmut Kaiser
//
//  Part of this code has been adopted from code published under the BSL by:
//
//  (C) Copyright 2005-7 Anthony Williams
//  (C) Copyright 2007 David Deakins
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_60DE7FDD_030E_4E7F_B47D_0C2AA53AE4BF)
#define OCTOPUS_60DE7FDD_030E_4E7F_B47D_0C2AA53AE4BF

#include <boost/atomic.hpp>
#include <boost/integer.hpp>
#include <boost/noncopyable.hpp>

namespace octopus
{

template <boost::uint64_t Bits>
struct atomic_bitset : boost::noncopyable
{
    typedef typename boost::uint_t<Bits>::exact storage_type;

    typedef boost::uint64_t size_type;

  private:
    boost::atomic<storage_type> bits_;

  public:
    atomic_bitset() : bits_(0) {}

    atomic_bitset& set()
    {
        bits_.store(~storage_type(0));
        return *this; 
    }

    atomic_bitset& set(size_type pos, bool val)
    {
        OCTOPUS_ASSERT(pos < Bits);

        if (!val)
            return reset(pos);

        storage_type const mask = 1 << pos;
        storage_type old = bits_.load(boost::memory_order_acquire);

        do {
            storage_type tmp = old;
            if (bits_.compare_exchange_strong(tmp, old | mask))
                break;
            old = tmp;
        } while (true);

        return *this; 
    }

    atomic_bitset& reset()
    {
        bits_.store(storage_type(0));
        return *this;
    }

    atomic_bitset& reset(size_type pos)
    {
        OCTOPUS_ASSERT(pos < Bits);

        storage_type const mask = 1 << pos;
        storage_type old = bits_.load(boost::memory_order_acquire);

        do {
            storage_type tmp = old;
            if (bits_.compare_exchange_strong(tmp, old & ~mask))
                break;
            old = tmp;
        } while (true);

        return *this; 
    }

    bool test(size_type pos)
    {
        OCTOPUS_ASSERT(pos < Bits);
        return bits_.load() & (1 << pos);
    }

    storage_type bits()
    {
        return bits_.load();
    }
}; 

}

#endif // OCTOPUS_60DE7FDD_030E_4E7F_B47D_0C2AA53AE4BF

