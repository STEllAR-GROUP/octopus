////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA)
#define OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA

#include <hpx/runtime/naming/name.hpp>
#include <hpx/lcos/future.hpp>

#include <octopus/child_index.hpp>

#include <boost/serialization/access.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT octree_client
{
  private:
    hpx::id_type gid_;

    friend struct octree_server;

    void create(hpx::id_type const& locality);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & gid_;
    }

  public:
    octree_client() : gid_(hpx::naming::invalid_id) {}

    octree_client(octree_client const& other) : gid_(other.gid_) {}

    octree_client(BOOST_RV_REF(octree_client) other) : gid_(other.gid_)
    {
        other.gid_.reset();
    }

    octree_client& operator=(BOOST_COPY_ASSIGN_REF(octree_client) other)
    {
        if (gid_ != other.gid_)
            gid_ = other.gid_;
        return *this;
    }

    octree_client& operator=(BOOST_RV_REF(octree_client) other)
    {
        if (gid_ != other.gid_)
        {
            gid_ = other.gid_;
            other.gid_.reset();
        }
        return *this;
    }

    octree_client& operator=(hpx::id_type const& gid)
    {
        if (gid_ != gid)
            gid_ = gid;
        return *this;
    }

    octree_client& operator=(BOOST_RV_REF(hpx::id_type) gid)
    {
        if (gid_ != gid)
            gid_ = gid;
        return *this;
    }

    hpx::id_type const& get_gid() const
    {
        return gid_;
    }

    ///////////////////////////////////////////////////////////////////////////
    void create_child(child_index idx);

    hpx::future<void> create_child_async(child_index idx);
};

}

#endif // OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA

