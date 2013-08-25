////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_1B486746_2FA0_455E_B9D5_84D457F3533C)
#define OCTOPUS_1B486746_2FA0_455E_B9D5_84D457F3533C

#include <octopus/assert.hpp>
#include <octopus/array.hpp>

#include <boost/move/move.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/access.hpp>

namespace octopus
{

template <
    typename T
  , boost::uint64_t TotalLength = OCTOPUS_GRID_NODE_LENGTH
                                * OCTOPUS_GRID_NODE_LENGTH
                                * OCTOPUS_GRID_NODE_LENGTH
    >
struct array3d
{
    typedef boost::uint64_t size_type;
    typedef boost::array<T, TotalLength> storage_type;

  private:
    size_type x_length_;
    size_type y_length_;
    size_type z_length_;
    boost::shared_ptr<storage_type> data_;

    BOOST_COPYABLE_AND_MOVABLE(array3d);

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
        return x + y * x_length_ + z * x_length_ * y_length_;
    }

  public:
    array3d() {}

    array3d(
        size_type length
      , T const& dflt = T()
        )
      : x_length_(length)
      , y_length_(length)
      , z_length_(length) 
      , data_(new storage_type)
    {
        OCTOPUS_ASSERT(    (x_length_ * y_length_ * z_length_)
                        == TotalLength);

        for (boost::uint64_t i = 0; i < data_->size(); ++i)
            (*data_)[i] = dflt;
    }

    array3d(
        size_type x_length
      , size_type y_length
      , size_type z_length
      , T const& dflt = T()
        )
      : x_length_(x_length)
      , y_length_(y_length)
      , z_length_(z_length) 
      , data_(new storage_type)
    {
        OCTOPUS_ASSERT(    (x_length_ * y_length_ * z_length_)
                        == TotalLength);

        for (boost::uint64_t i = 0; i < data_->size(); ++i)
            (*data_)[i] = dflt;
    }

    array3d(array3d const& other)
      : x_length_(other.x_length_)
      , y_length_(other.y_length_)
      , z_length_(other.z_length_) 
      , data_(new storage_type)
    {
        for (boost::uint64_t i = 0; i < data_->size(); ++i)
            (*data_)[i] = (*other.data_)[i];
    }

    array3d(BOOST_RV_REF(array3d) other)
      : x_length_(other.x_length_)
      , y_length_(other.y_length_)
      , z_length_(other.z_length_) 
      , data_(boost::move(other.data_))
    {}

    array3d& operator=(BOOST_COPY_ASSIGN_REF(array3d) other)
    {
        x_length_ = other.x_length_;
        y_length_ = other.y_length_;
        z_length_ = other.z_length_;

        if (!data_)
            data_.reset(new storage_type);

        for (boost::uint64_t i = 0; i < data_->size(); ++i)
            (*data_)[i] = (*other.data_)[i];

        return *this;
    }

    array3d& operator=(BOOST_RV_REF(array3d) other)
    {
        x_length_ = other.x_length_;
        y_length_ = other.y_length_;
        z_length_ = other.z_length_;
        data_ = boost::move(other.data_);

        return *this;
    }

    array3d& operator=(T const& value)
    {
        if (!data_)
            data_.reset(new storage_type);

        for (size_type i = 0; i < data_->size(); ++i)
            (*data_)[i] = value;

        return *this;
    }

    size_type size() const
    {
        return TotalLength;
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

    T&
    operator()(size_type x, size_type y, size_type z)
    {
        return (*data_)[index(x, y, z)];
    }

    T const&
    operator()(size_type x, size_type y, size_type z) const
    {
        return (*data_)[index(x, y, z)];
    }

    T&
    operator()(array<size_type, 3> const& xyz)
    {
        return (*this)(xyz[0], xyz[1], xyz[2]);
    }

    T const&
    operator()(array<size_type, 3> const& xyz) const
    {
        return (*this)(xyz[0], xyz[1], xyz[2]);
    }
};

}

#endif // OCTOPUS_1B486746_2FA0_455E_B9D5_84D457F3533C

