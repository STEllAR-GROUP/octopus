////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_0CBAD952_ED85_4D22_84DF_6AB16A34E518)
#define OCTOPUS_0CBAD952_ED85_4D22_84DF_6AB16A34E518

#include <octopus/assert.hpp>
#include <octopus/face.hpp>

#include <boost/serialization/access.hpp>
#include <boost/cstdint.hpp>
#include <boost/array.hpp>

#include <ostream>

namespace octopus
{

struct child_index
{
    child_index() : packed_(0) {}

    child_index(boost::uint64_t packed)
      : packed_(packed)
    {
        OCTOPUS_ASSERT_MSG(7 >= packed, "packed index is out of range"); 
    }

    child_index(boost::uint64_t x, boost::uint64_t y, boost::uint64_t z)
      : packed_(x | (y << 1) | (z << 2))
    {
        OCTOPUS_ASSERT_MSG(1 >= x, "x-index is out of range"); 
        OCTOPUS_ASSERT_MSG(1 >= y, "y-index is out of range"); 
        OCTOPUS_ASSERT_MSG(1 >= z, "z-index is out of range"); 
    }

    child_index(child_index const& other)
      : packed_(other.packed_)
    {}

    child_index operator=(child_index const& other)
    {
        packed_ = other.packed_;
        return *this;
    }

    operator boost::uint64_t() const
    {
        return packed_;
    }

    bool x() const
    {
        return (packed_ >> 0) & 1;
    }

    bool y() const
    {
        return (packed_ >> 1) & 1;
    }

    bool z() const
    {
        return (packed_ >> 2) & 1;
    }

    void set_x(boost::uint64_t x_)
    {
        OCTOPUS_ASSERT_MSG(1 >= x_, "x-index is out of range"); 
        if (0 == x_)
            packed_ &= ~1;
        else
            packed_ |= 1;
    }

    void set_y(boost::uint64_t y_)
    {
        OCTOPUS_ASSERT_MSG(1 >= y_, "y-indey is out of range"); 
        if (0 == y_)
            packed_ &= ~2;
        else
            packed_ |= 2;
    }

    void set_z(boost::uint64_t z_)
    {
        OCTOPUS_ASSERT_MSG(1 >= z_, "z-indez is out of range"); 
        if (0 == z_)
            packed_ &= ~4;
        else
            packed_ |= 4;
    }

    boost::array<boost::int64_t, 3> array() const
    {
        boost::array<boost::int64_t, 3> a;
        a[0] = x();
        a[1] = y();
        a[2] = z();
        return a; 
    }  

  private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & packed_;
    }

    boost::uint8_t packed_;
};

inline std::ostream& operator<<(std::ostream& os, child_index idx)
{
    os << "(" << idx.x() << ", " << idx.y() << ", " << idx.z() << ")";
    return os;
}

OCTOPUS_EXPORT child_index invert(face f, child_index const& idx);

}

#endif // OCTOPUS_0CBAD952_ED85_4D22_84DF_6AB16A34E518

