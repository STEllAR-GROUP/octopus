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

#include <boost/serialization/access.hpp>
#include <boost/cstdint.hpp>

#include <ostream>

namespace octopus
{

struct child_index
{
    child_index() : packed_() {}

    child_index(std::size_t x_, std::size_t y_, std::size_t z_)
      : packed_(x_ | (y_ << 1) | (z_ << 2))
    {
        OCTOPUS_ASSERT_MSG(1 >= x_, "x-index is out of range"); 
        OCTOPUS_ASSERT_MSG(1 >= y_, "y-index is out of range"); 
        OCTOPUS_ASSERT_MSG(1 >= z_, "z-index is out of range"); 
    }

    child_index(child_index const& other)
    {
        *this = other;
    }

    child_index operator=(child_index const& other)
    {
        packed_ = other.packed_;
        return *this;
    }

    operator std::size_t() const
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
        return (packed_ >> 2) & 0x1;
    }

    void set_x(std::size_t x_)
    {
        OCTOPUS_ASSERT_MSG(1 >= x_, "x-index is out of range"); 
        if (0 == x_)
            packed_ &= ~1;
        else
            packed_ |= 1;
    }

    void set_y(std::size_t y_)
    {
        OCTOPUS_ASSERT_MSG(1 >= y_, "y-indey is out of range"); 
        if (0 == y_)
            packed_ &= ~2;
        else
            packed_ |= 2;
    }

    void set_z(std::size_t z_)
    {
        OCTOPUS_ASSERT_MSG(1 >= z_, "z-indez is out of range"); 
        if (0 == z_)
            packed_ &= ~4;
        else
            packed_ |= 4;
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

}

#endif // OCTOPUS_0CBAD952_ED85_4D22_84DF_6AB16A34E518

