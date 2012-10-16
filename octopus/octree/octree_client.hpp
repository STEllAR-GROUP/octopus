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

#include <octopus/octree/octree_init_data.hpp>
#include <octopus/child_index.hpp>
#include <octopus/face.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT octree_server;

/// The set of types in our type-punning system. We call these types "kinds",
/// to distinguish them from C++ types.
// REVIEW: Should this have a serialization version?
enum boundary_kind
{
    invalid_boundary
  , real_boundary
  , physical_boundary
  , amr_boundary 
};

inline std::ostream& operator<<(std::ostream& os, boundary_kind k)
{
    switch (k)
    {
        case real_boundary:     os << "real_boundary"; break; 
        case physical_boundary: os << "physical_boundary"; break; 
        case amr_boundary:      os << "amr_boundary"; break; 
        default: os << "invalid_boundary"; break; 
    }
    return os;
}

}

///////////////////////////////////////////////////////////////////////////////
namespace boost { namespace serialization
{
    template <typename Archive>
    void save(Archive& ar, octopus::boundary_kind const& k, const unsigned int)
    {
        boost::uint8_t tmp(k);
        ar & tmp; 
    }

    template <typename Archive>
    void load(Archive& ar, octopus::boundary_kind& k, const unsigned int)
    {
        boost::uint8_t tmp;
        ar & tmp; 
        OCTOPUS_ASSERT_FMT_MSG(octopus::invalid_boundary > k,
                               "invalid face deserialized, face(%1%)",
                               boost::uint16_t(tmp));
        k = tmp; 
    }
}}

BOOST_SERIALIZATION_SPLIT_FREE(octopus::boundary_kind);

namespace octopus
{

struct OCTOPUS_EXPORT octree_client
{
  private:
    hpx::id_type gid_;
    boundary_kind kind_;

    // FIXME: This is only used by non-real boundaries.
    face face_;

    // FIXME: This is only used for physical boundaries, optimize.
    vector3d<std::vector<double> > map_;
    bool reflect_;
    boost::uint64_t direction_; ///< 0 == x, 1 == y, 2 == z 

    // FIXME: This is only used for AMR boundaries, optimize.
    boost::array<boost::int64_t, 3> offset_; ///< Relative offset.
    
    BOOST_COPYABLE_AND_MOVABLE(octree_client);

    friend struct octree_server;

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & gid_;
        ar & kind_;
        if (!real())
        {
            ar & face_;
            ar & map_;
            ar & reflect_;
            ar & direction_;
            ar & offset_;
        }
    }

    /// \brief Create a client for a real grid node.
    ///
    /// \note Internal use only; GIDs are not exposed to the user.
    octree_client(hpx::id_type const& gid)
      : gid_(gid)
      , kind_(real_boundary)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {}

    /// \brief Create a client for a real grid node.
    ///
    /// \note Internal use only; GIDs are not exposed to the user.
    octree_client(BOOST_RV_REF(hpx::id_type) gid)
      : gid_(gid)
      , kind_(real_boundary)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {}

    octree_client(hpx::id_type const& parent, boundary_kind kind)
      : gid_(parent)
      , kind_(kind)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {
        OCTOPUS_ASSERT(kind != real_boundary);
    }

    octree_client(BOOST_RV_REF(hpx::id_type) parent, boundary_kind kind)
      : gid_(parent)
      , kind_(kind)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {
        OCTOPUS_ASSERT(kind != real_boundary);
    }

    /// \brief Assign the GID of a real grid node to this client. 
    ///
    /// \note Internal use only; GIDs are not exposed to the user.
    octree_client& operator=(hpx::id_type const& gid)
    {
        gid_ = gid;
        kind_ = real_boundary;
        face_ = face();
        map_.clear();
        reflect_ = false;
        direction_ = 0;
        offset_[0] = 0;
        offset_[1] = 0;
        offset_[2] = 0;
        return *this;
    }

    /// \brief Assign the GID of a real grid node to this client. 
    ///
    /// \note Internal use only; GIDs are not exposed to the user.
    octree_client& operator=(BOOST_RV_REF(hpx::id_type) gid)
    {
        gid_ = gid;
        kind_ = real_boundary;
        face_ = face();
        map_.clear();
        reflect_ = false;
        direction_ = 0;
        offset_[0] = 0;
        offset_[1] = 0;
        offset_[2] = 0;
        return *this;
    }

    /// \brief Get the GID this client is holding. 
    ///
    /// \note Internal use only; GIDs are not exposed to the user.
    hpx::id_type const& get_gid() const
    {
        return gid_;
    }

    void ensure_real()
    {
        OCTOPUS_ASSERT_FMT_MSG(kind_ == real_boundary,
            "illegal operation for %1% client, expected real boundary",
            kind_);
    }

  public:
    octree_client()
      : gid_(hpx::naming::invalid_id)
      , kind_(invalid_boundary)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {}

    octree_client(octree_client const& other)
      : gid_(other.gid_)
      , kind_(other.kind_)
      , face_(other.face_)
      , map_(other.map_)
      , reflect_(other.reflect_)
      , direction_(other.direction_)
      , offset_(other.offset_)
    {}

    octree_client(BOOST_RV_REF(octree_client) other)
      : gid_(other.gid_)
      , kind_(other.kind_)
      , face_(other.face_)
      , map_(other.map_)
      , reflect_(other.reflect_)
      , direction_(other.direction_)
      , offset_(other.offset_)
    {}

    octree_client(octree_client const& parent, boundary_kind kind)
      : gid_(parent.gid_)
      , kind_(kind)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {
        OCTOPUS_ASSERT(kind != real_boundary);
    }

    octree_client(BOOST_RV_REF(octree_client) parent, boundary_kind kind)
      : gid_(parent.gid_)
      , kind_(kind)
      , face_()
      , map_()
      , reflect_()
      , direction_()
      , offset_()
    {
        OCTOPUS_ASSERT(kind != real_boundary);
    }

    octree_client& operator=(BOOST_COPY_ASSIGN_REF(octree_client) other)
    {
        gid_ = other.gid_;
        kind_ = other.kind_;        
        face_ = other.face_;
        map_ = other.map_;
        reflect_ = other.reflect_;
        direction_ = other.direction_;
        offset_ = other.offset_;
        return *this;
    }

    octree_client& operator=(BOOST_RV_REF(octree_client) other)
    {
        gid_ = other.gid_;
        kind_ = other.kind_;        
        face_ = other.face_;
        map_ = other.map_;
        reflect_ = other.reflect_;
        direction_ = other.direction_;
        offset_ = other.offset_;
        return *this;
    }

    operator hpx::util::safe_bool<octree_client>::result_type() const
    {
        return hpx::util::safe_bool<octree_client>()(gid_);
    }

    friend bool operator==(octree_client const& lhs, octree_client const& rhs) 
    {
        return lhs.gid_ == rhs.gid_
            && lhs.kind_ == rhs.kind_; 
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
        return !(lhs == rhs);
    }

    friend bool operator!=(octree_client const& lhs, hpx::id_type const& rhs) 
    {
        return !(lhs == rhs);
    }

    friend bool operator!=(hpx::id_type const& lhs, octree_client const& rhs) 
    {
        return !(lhs == rhs);
    }
    
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Returns true if this client represents a real octree (as opposed
    ///        an AMR or physical boundary).
    bool real() const 
    {
        switch (kind_)
        {
            case real_boundary:
                return true;
            case physical_boundary:
            case amr_boundary:
                return false; 
            default:
                break;
        };

        OCTOPUS_ASSERT_MSG(false, "invalid boundary kind");
        return false;
    }

    /// \brief Return the boundary_kind of this client.
    boundary_kind kind() const
    {
        return kind_;
    }

    ///////////////////////////////////////////////////////////////////////////
    // {{{ create_root
    void create_root(
        hpx::id_type const& locality
      , octree_init_data const& init   
        );

    void create_root(
        hpx::id_type const& locality
      , BOOST_RV_REF(octree_init_data) init   
        );
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ create_child
    void create_child(
        child_index kid
        )
    {
        create_child_async(kid).get(); 
    }

    hpx::future<void> create_child_async(
        child_index kid
        );
    // }}}

  private:
    void set_sibling_for_amr_boundary(
        face f
      , octree_client const& sib 
        );

    void tie_sibling_for_amr_boundary(
        face f
      , octree_client const& sib 
        );

    void set_sibling_for_physical_boundary(
        face f
      , octree_client const& sib 
        );

    void tie_sibling_for_physical_boundary(
        face f
      , octree_client const& sib 
        );
    
  public:
    ///////////////////////////////////////////////////////////////////////////
    // {{{ set_sibling
    void set_sibling(
        face f
      , octree_client const& sib 
        );

    void set_sibling_push(
        face f
      , octree_client const& sib 
        );
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ tie_sibling
    void tie_sibling(
        face f
      , octree_client const& sib 
        );

    void tie_sibling_push(
        face f
      , octree_client const& sib 
        );
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ tie_child_sibling
    void set_child_sibling(
        child_index kid
      , face f
      , octree_client const& sib 
        );

    void set_child_sibling_push(
        child_index kid
      , face f
      , octree_client const& sib 
        );
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ tie_child_sibling
    void tie_child_sibling(
        child_index kid
      , face f
      , octree_client const& sib 
        );

    void tie_child_sibling_push(
        child_index kid
      , face f
      , octree_client const& sib 
        );
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ get_siblings
    boost::array<octree_client, 6> get_siblings()
    {
        return get_siblings_async().get();
    }

    hpx::future<boost::array<octree_client, 6> > get_siblings_async();
    // }}}

    ///////////////////////////////////////////////////////////////////////////
  private:
    vector3d<std::vector<double> > interpolate(
        face f
        );

    vector3d<std::vector<double> > mirror_or_outflow(
        face f
        );

  public:
    // {{{ receive_ghost_zones
    void receive_ghost_zones()
    {
        receive_ghost_zones_async().get();
    }

    hpx::future<void> receive_ghost_zones_async();
    // }}}

    // {{{ send_ghost_zone
    vector3d<std::vector<double> > send_ghost_zone(
        face f
        );

    hpx::future<vector3d<std::vector<double> > > send_ghost_zone_async(
        face f
        );
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ apply
    void apply(
        hpx::util::function<void(octree_server&)> const& f
      , boost::uint64_t minimum_level
        )
    {
        return apply_async(f, minimum_level).get();
    }

    hpx::future<void> apply_async(
        hpx::util::function<void(octree_server&)> const& f
      , boost::uint64_t minimum_level
        );

    void apply_push(
        hpx::util::function<void(octree_server&)> const& f
      , boost::uint64_t minimum_level
        );
    // }}} 

    // NOTE: (to self) Keep the order the same as octree_server please.
};

}

#endif // OCTOPUS_17096997_37B3_4F9E_80F3_4C964006BFAA

