////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA)
#define OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA

#include <hpx/runtime/naming/name.hpp>
#include <hpx/lcos/future.hpp>

#include <octopus/array1d.hpp>
#include <octopus/child_index.hpp>
#include <octopus/face.hpp>

#include <boost/serialization/access.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT octree_client
{
  private:
    hpx::id_type gid_;

    friend struct octree_server;

    void create(
        hpx::id_type const& locality
      , boost::uint64_t level
      , array1d<boost::uint64_t, 3> const& location
        );

    // NOTE: Does not set the GID of this client.
    hpx::future<hpx::id_type, hpx::naming::gid_type> create_async(
        hpx::id_type const& locality
      , boost::uint64_t level
      , array1d<boost::uint64_t, 3> const& location
        ) const;

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & gid_;
    }

  public:
    octree_client() : gid_(hpx::naming::invalid_id) {}

    octree_client(octree_client const& other) : gid_(other.gid_) {}

    octree_client(hpx::id_type const& gid) : gid_(gid) {}

    octree_client(BOOST_RV_REF(octree_client) other) : gid_(other.gid_)
    {
        other.gid_ = hpx::naming::invalid_id;
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
            other.gid_ = hpx::naming::invalid_id;
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

    operator hpx::util::safe_bool<octree_client>::result_type() const
    {
        return hpx::util::safe_bool<octree_client>()(gid_);
    }

    friend bool operator==(octree_client const& lhs, octree_client const& rhs) 
    {
        return lhs.gid_ == rhs.gid_;
    }

    friend bool operator==(octree_client const& lhs, hpx::id_type const& rhs) 
    {
        return lhs.gid_ == rhs;
    }

    friend bool operator==(hpx::id_type const& lhs, octree_client const& rhs) 
    {
        return lhs == rhs.gid_;
    }

    friend bool operator!=(octree_client const& lhs, octree_client const& rhs) 
    {
        return lhs.gid_ != rhs.gid_;
    }

    friend bool operator!=(octree_client const& lhs, hpx::id_type const& rhs) 
    {
        return lhs.gid_ != rhs;
    }

    friend bool operator!=(hpx::id_type const& lhs, octree_client const& rhs) 
    {
        return lhs != rhs.gid_;
    }

    void create_root(
        hpx::id_type const& gid
        );

    ///////////////////////////////////////////////////////////////////////////
    void create_child(
        child_index kid
        )
    {
        create_child_async(kid).get(); 
    }

    hpx::future<void> create_child_async(
        child_index kid
        );

    ///////////////////////////////////////////////////////////////////////////
  private:
    void set_sibling(
        boost::uint8_t f
      , octree_client const& sib 
        );

    void set_sibling_push(
        boost::uint8_t f
      , octree_client const& sib 
        );

  public:
    void set_sibling(
        face f
      , octree_client const& sib 
        )
    {
        set_sibling(boost::uint8_t(f), sib);
    }

    void set_sibling_push(
        face f
      , octree_client const& sib 
        )
    {
        return set_sibling_push(boost::uint8_t(f), sib);
    }

    ///////////////////////////////////////////////////////////////////////////
  private:
    void tie_sibling(
        boost::uint8_t f
      , octree_client const& sib 
        );

    void tie_sibling_push(
        boost::uint8_t f
      , octree_client const& sib 
        );

  public:
    void tie_sibling(
        face f
      , octree_client const& sib 
        )
    {
        set_sibling(boost::uint8_t(f), sib);
    }

    void tie_sibling_push(
        face f
      , octree_client const& sib 
        )
    {
        set_sibling_push(boost::uint8_t(f), sib);
    }

    ///////////////////////////////////////////////////////////////////////////
  private:
    void set_child_sibling(
        child_index kid
      , boost::uint8_t f
      , octree_client const& sib 
        );

    void set_child_sibling_push(
        child_index kid
      , boost::uint8_t f
      , octree_client const& sib 
        );

  public:
    void set_child_sibling(
        child_index kid
      , face f
      , octree_client const& sib 
        )
    {
        set_child_sibling(kid, boost::uint8_t(f), sib);
    }

    void set_child_sibling_push(
        child_index kid
      , face f
      , octree_client const& sib 
        )
    {
        set_child_sibling_push(kid, boost::uint8_t(f), sib);
    }

    ///////////////////////////////////////////////////////////////////////////
  private:
    void tie_child_sibling(
        child_index kid
      , boost::uint8_t f
      , octree_client const& sib 
        );

    void tie_child_sibling_push(
        child_index kid
      , boost::uint8_t f
      , octree_client const& sib 
        );

  public:
    void tie_child_sibling(
        child_index kid
      , face f
      , octree_client const& sib 
        )
    {
        tie_child_sibling(kid, boost::uint8_t(f), sib);
    }

    void tie_child_sibling_push(
        child_index kid
      , face f
      , octree_client const& sib 
        )
    {
        tie_child_sibling_push(kid, boost::uint8_t(f), sib);
    }
};

}

#endif // OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA

