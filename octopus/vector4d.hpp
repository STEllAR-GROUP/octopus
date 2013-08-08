////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_1D91DF9D_CF4E_4616_8E65_AE959E8533FD)
#define OCTOPUS_1D91DF9D_CF4E_4616_8E65_AE959E8533FD

#include <octopus/assert.hpp>
#include <octopus/array.hpp>

#include <boost/move/move.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/access.hpp>

#include <vector>

namespace octopus
{

template <typename T, boost::uint64_t SLength = OCTOPUS_STATE_SIZE>
struct vector4d
{
    typedef boost::uint64_t size_type;

  private:
    size_type x_length_;
    size_type y_length_;
    size_type z_length_;
    std::vector<T> data_;

    BOOST_COPYABLE_AND_MOVABLE(vector4d);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & x_length_ & y_length_ & z_length_ & data_;
    }

    size_type index(size_type x, size_type y, size_type z) const
    {
        OCTOPUS_ASSERT_FMT_MSG(x < x_length_,
            "x coordinate (%1%) is larger than the x length (%2%)",
            x % x_length_);  
        OCTOPUS_ASSERT_FMT_MSG(y < y_length_,
            "y coordinate (%1%) is larger than the y length (%2%)",
            y % y_length_);  
        OCTOPUS_ASSERT_FMT_MSG(z < z_length_,
            "z coordinate (%1%) is larger than the z length (%2%)",
            z % z_length_);  
        return x * SLength
             + y * x_length_ * SLength
             + z * x_length_ * y_length_ * SLength;
    }

  public:
    vector4d() {}

    vector4d(
        size_type length
      , T const& dflt = T()
        )
      : x_length_(length)
      , y_length_(length)
      , z_length_(length) 
      , data_(x_length_ * y_length_ * z_length_ * SLength, dflt)
    {}

    vector4d(
        size_type x_length
      , size_type y_length
      , size_type z_length
      , T const& dflt = T()
        )
      : x_length_(x_length)
      , y_length_(y_length)
      , z_length_(z_length) 
      , data_(x_length_ * y_length_ * z_length_ * SLength, dflt)
    {}

    vector4d(vector4d const& other)
      : x_length_(other.x_length_)
      , y_length_(other.y_length_)
      , z_length_(other.z_length_) 
      , data_(other.data_)
    {}

    vector4d(BOOST_RV_REF(vector4d) other)
      : x_length_(other.x_length_)
      , y_length_(other.y_length_)
      , z_length_(other.z_length_) 
      , data_(boost::move(other.data_))
    {
        other.clear();
    }

    vector4d& operator=(BOOST_COPY_ASSIGN_REF(vector4d) other)
    {
        x_length_ = other.x_length_;
        y_length_ = other.y_length_;
        z_length_ = other.z_length_;
        data_ = other.data_;
        return *this;
    }

    vector4d& operator=(BOOST_RV_REF(vector4d) other)
    {
        x_length_ = other.x_length_;
        y_length_ = other.y_length_;
        z_length_ = other.z_length_;
        data_ = boost::move(other.data_);

        other.clear();

        return *this;
    }

    vector4d& operator=(T const& value)
    {
        for (size_type i = 0; i < data_.size(); ++i)
            data_[i] = value;
        return *this;
    }

    size_type size() const
    {
        return data_.size() / SLength;
    } 

    size_type x_length() const
    {
        return x_length_; 
    } 

    size_type y_length() const
    {
        return y_length_; 
    } 

    size_type z_length() const
    {
        return z_length_; 
    } 

    size_type s_length() const
    {
        return SLength; 
    } 

    void resize(
        size_type length
      , T const& dflt = T()
        )
    {
        x_length_ = length;
        y_length_ = length;
        z_length_ = length;
        data_.resize(x_length_ * y_length_ * z_length_ * SLength, dflt);    
    }

    void resize(
        size_type x_length
      , size_type y_length
      , size_type z_length
      , T const& dflt = T()
        )
    {
        x_length_ = x_length;
        y_length_ = y_length;
        z_length_ = z_length;
        data_.resize(x_length_ * y_length_ * z_length_ * SLength, dflt);    
    }

    void clear()
    {
        x_length_ = 0;
        y_length_ = 0;
        z_length_ = 0;
        data_.clear();
    }

    array<T, SLength>&
    operator()(size_type x, size_type y, size_type z)
    {
        return *((array<T, SLength>*) &data_[index(x, y, z)]);
    }

    array<T, SLength> const&
    operator()(size_type x, size_type y, size_type z) const
    {
        return *((array<T, SLength>*) &data_[index(x, y, z)]);
    }

    array<T, SLength>&
    operator()(array<size_type, 3> const& xyz)
    {
        return (*this)(xyz[0], xyz[1], xyz[2]);
    }

    array<T, SLength> const&
    operator()(array<size_type, 3> const& xyz) const
    {
        return (*this)(xyz[0], xyz[1], xyz[2]);
    }
};

}

#endif // OCTOPUS_820FC460_190E_4E57_8CB7_CB0880FC3E28

