////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_820FC460_190E_4E57_8CB7_CB0880FC3E28)
#define OCTOPUS_820FC460_190E_4E57_8CB7_CB0880FC3E28

#include <octopus/assert.hpp>

#include <boost/move/move.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/access.hpp>

#include <vector>

// TODO: I have limited this to a cube currently because that's all we need. 

namespace octopus
{

template <typename T>
struct vector3d
{
  private:
    std::size_t dimension_;
    std::vector<T> data_;

    BOOST_COPYABLE_AND_MOVABLE(vector3d);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & dimension_ & data_;
    }

  public:
    vector3d() {}

    vector3d(std::size_t dimension, T dflt = T())
      : dimension_(dimension), data_(dimension_ * dimension_ * dimension_, dflt)
    {}

    vector3d(vector3d const& other)
      : dimension_(other.dimension_), data_(other.data_)
    {}

    vector3d(BOOST_RV_REF(vector3d) other)
      : dimension_(other.dimension_), data_(other.data_)
    {}

    vector3d& operator=(BOOST_COPY_ASSIGN_REF(vector3d) other)
    {
        dimension_ = other.dimension_;
        data_ = other.data_;
        return *this;
    }

    vector3d& operator=(BOOST_RV_REF(vector3d) other)
    {
        dimension_ = other.dimension_;
        data_ = other.data_;
        return *this;
    }

    vector3d& operator=(T const& value)
    {
        for (std::size_t i = 0; i < data_.size(); ++i)
            data_[i] = value;
        return *this;
    }

    T& operator()(std::size_t x, std::size_t y, std::size_t z)
    {
        OCTOPUS_ASSERT_FMT_MSG(x < dimension_,
            "x coordinate (%1%) is smaller than the x dimension (%2%)",
            x % dimension_);  
        OCTOPUS_ASSERT_FMT_MSG(y < dimension_,
            "y coordinate (%1%) is smaller than the y dimension (%2%)",
            y % dimension_);  
        OCTOPUS_ASSERT_FMT_MSG(z < dimension_,
            "z coordinate (%1%) is smaller than the z dimension (%2%)",
            z % dimension_);  
        return data_[x + y * dimension_ + z * dimension_ * dimension_];
    }

    T operator()(std::size_t x, std::size_t y, std::size_t z) const
    {
        OCTOPUS_ASSERT_FMT_MSG(x < dimension_,
            "x coordinate (%1%) is smaller than the x dimension (%2%)",
            x % dimension_);  
        OCTOPUS_ASSERT_FMT_MSG(y < dimension_,
            "y coordinate (%1%) is smaller than the y dimension (%2%)",
            y % dimension_);  
        OCTOPUS_ASSERT_FMT_MSG(z < dimension_,
            "z coordinate (%1%) is smaller than the z dimension (%2%)",
            z % dimension_);  
        return data_[x + y * dimension_ + z * dimension_ * dimension_];
    }

    vector3d& operator+=(T const& v)
    {
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] += v;
        return *this;
    }

    vector3d& operator-=(T const& v)
    {
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] -= v;
        return *this;
    }

    vector3d& operator*=(T const& v)
    {
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] *= v;
        return *this;
    }

    vector3d& operator/=(T const& v)
    {
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] /= v;
        return *this;
    }

    vector3d& operator+=(vector3d const& v)
    { 
        OCTOPUS_ASSERT_MSG(v.dimension_ != dimension_,
            "dimensions do not match");
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] += v[i];
        return *this;
    }

    vector3d& operator-=(vector3d const& v)
    {
        OCTOPUS_ASSERT_MSG(v.dimension_ != dimension_,
            "dimensions do not match");
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] -= v[i];
        return *this;
    }

    vector3d& operator*=(vector3d const& v)
    {
        OCTOPUS_ASSERT_MSG(v.dimension_ != dimension_,
            "dimensions do not match");
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] *= v[i];
        return *this;
    }

    vector3d& operator/=(vector3d const& v)
    {
        OCTOPUS_ASSERT_MSG(v.dimension_ != dimension_,
            "dimensions do not match");
        for (std::size_t i = 0; i < data_.size(); ++i) 
            data_[i] /= v[i];
        return *this;
    }

    vector3d operator+(vector3d const& v) const
    {
        return (vector3d(*this) += v);
    }

    vector3d operator-(vector3d const& v) const
    {
        return (vector3d(*this) -= v);
    }

    vector3d operator*(vector3d const& v) const
    {
        return (vector3d(*this) *= v);
    }

    vector3d operator/(vector3d const& v) const
    {
        return (vector3d(*this) /= v);
    }

    vector3d operator+(T const& v) const
    {
        return (vector3d(*this) += v);
    }

    vector3d operator-(T const& v) const
    {
        return (vector3d(*this) -= v);
    }

    vector3d operator*(T const& v) const
    {
        return (vector3d(*this) *= v);
    }

    vector3d operator/(T const& v) const
    {
        return (vector3d(*this) /= v);
    }
};

}

#endif // OCTOPUS_820FC460_190E_4E57_8CB7_CB0880FC3E28

