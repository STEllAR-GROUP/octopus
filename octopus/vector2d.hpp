////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_6EEC1D66_A902_4089_9CD9_54B9F4FCE90B)
#define OCTOPUS_6EEC1D66_A902_4089_9CD9_54B9F4FCE90B

#include <octopus/assert.hpp>
#include <octopus/array.hpp>

#include <boost/move/move.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/access.hpp>

#include <vector>

namespace octopus
{

template <typename T, boost::uint64_t SLength = OCTOPUS_STATE_SIZE>
struct vector2d
{
    typedef boost::uint64_t size_type;

  private:
    size_type length_;
    std::vector<T> data_;

    BOOST_COPYABLE_AND_MOVABLE(vector2d);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & length_ & data_;
    }

    size_type index(size_type x) const
    {
        OCTOPUS_ASSERT_FMT_MSG(x < length_,
            "coordinate (%1%) is larger than the length (%2%)",
            x % length_);  
        return x * SLength;
    }

  public:
    vector2d() {}

    vector2d(
        size_type length
      , T const& dflt = T()
        )
      : length_(length)
      , data_(length_ * SLength, dflt)
    {}

    vector2d(vector2d const& other)
      : length_(other.length_)
      , data_(other.data_)
    {}

    vector2d(BOOST_RV_REF(vector2d) other)
      : length_(other.length_)
      , data_(boost::move(other.data_))
    {
        //other.clear();
    }

    vector2d& operator=(BOOST_COPY_ASSIGN_REF(vector2d) other)
    {
        length_ = other.length_;
        data_ = other.data_;
        return *this;
    }

    vector2d& operator=(BOOST_RV_REF(vector2d) other)
    {
        length_ = other.length_;
        data_ = boost::move(other.data_);

        //other.clear();

        return *this;
    }

    vector2d& operator=(T const& value)
    {
        for (size_type i = 0; i < data_.size(); ++i)
            data_[i] = value;
        return *this;
    }

    size_type size() const
    {
        return data_.size() / SLength;
    } 

    void resize(
        size_type length
      , T const& dflt = T()
        )
    {
        length_ = length;
        data_.resize(length_ * SLength, dflt);    
    }

    void clear()
    {
        length_ = 0;
        data_.clear();
    }

    array<T, SLength>&
    operator()(size_type x)
    {
        return *((array<T, SLength>*) &data_[index(x)]);
    }

    array<T, SLength> const&
    operator()(size_type x) const
    {
        return *((array<T, SLength>*) &data_[index(x)]);
    }

    array<T, SLength>&
    operator[](size_type x)
    {
        return *((array<T, SLength>*) &data_[index(x)]);
    }

    array<T, SLength> const&
    operator[](size_type x) const
    {
        return *((array<T, SLength>*) &data_[index(x)]);
    }
};

}

#endif // OCTOPUS_6EEC1D66_A902_4089_9CD9_54B9F4FCE90B

