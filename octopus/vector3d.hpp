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
#include <octopus/operators/std_vector_arithmetic.hpp>

#include <boost/move/move.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/access.hpp>

#include <vector>

namespace octopus
{

template <typename T>
struct vector3d
{
    typedef boost::uint64_t size_type;

  private:
    size_type x_length_;
    size_type y_length_;
    size_type z_length_;
    std::vector<T> data_;

    BOOST_COPYABLE_AND_MOVABLE(vector3d);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & x_length_ & y_length_ & z_length_ & data_;
    }

  public:
    vector3d() {}

    vector3d(
        size_type length
      , T dflt = T()
        )
      : x_length_(length)
      , y_length_(length)
      , z_length_(length) 
      , data_(x_length_ * y_length_ * z_length_, dflt)
    {}

    vector3d(
        size_type x_length
      , size_type y_length
      , size_type z_length
      , T dflt = T()
        )
      : x_length_(x_length)
      , y_length_(y_length)
      , z_length_(z_length) 
      , data_(x_length_ * y_length_ * z_length_, dflt)
    {}

    vector3d(vector3d const& other)
      : x_length_(other.x_length_)
      , y_length_(other.y_length_)
      , z_length_(other.z_length_) 
      , data_(other.data_)
    {}

    vector3d(BOOST_RV_REF(vector3d) other)
      : x_length_(other.x_length_)
      , y_length_(other.y_length_)
      , z_length_(other.z_length_) 
      , data_(other.data_)
    {}

    vector3d& operator=(BOOST_COPY_ASSIGN_REF(vector3d) other)
    {
        x_length_ = other.x_length_;
        y_length_ = other.y_length_;
        z_length_ = other.z_length_;
        data_ = other.data_;
        return *this;
    }

    vector3d& operator=(BOOST_RV_REF(vector3d) other)
    {
        x_length_ = other.x_length_;
        y_length_ = other.y_length_;
        z_length_ = other.z_length_;
        data_ = other.data_;
        return *this;
    }

    vector3d& operator=(T const& value)
    {
        for (size_type i = 0; i < data_.size(); ++i)
            data_[i] = value;
        return *this;
    }

    size_type size() const
    {
        return data_.size();
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

    void resize(
        size_type length
      , T dflt = T()
        )
    {
        x_length_ = length;
        y_length_ = length;
        z_length_ = length;
        data_.resize(x_length_ * y_length_ * z_length_, dflt);    
    }

    void resize(
        size_type x_length
      , size_type y_length
      , size_type z_length
      , T dflt = T()
        )
    {
        x_length_ = x_length;
        y_length_ = y_length;
        z_length_ = z_length;
        data_.resize(x_length * y_length * z_length, dflt);    
    }

    void clear()
    {
        x_length_ = 0;
        y_length_ = 0;
        z_length_ = 0;
        data_.clear();
    }

    T& operator()(size_type x, size_type y, size_type z)
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
        return data_[x + y * x_length_ + z * x_length_ * y_length_];
    }

    T const& operator()(size_type x, size_type y, size_type z) const
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
        return data_[x + y * x_length_ + z * x_length_ * y_length_];
    }

    T& operator()(boost::array<size_type, 3> const& xyz)
    {
        return (*this)(xyz[0], xyz[1], xyz[2]);
    }

    T const& operator()(boost::array<size_type, 3> const& xyz) const
    {
        return (*this)(xyz[0], xyz[1], xyz[2]);
    }

    T& operator[](size_type idx)
    {
        OCTOPUS_ASSERT_FMT_MSG(idx < data_.size(),
            "index (%1%) is larger than the data size (%2%)",
            idx % data_.size()); 
        return data_[idx];
    }

    T const& operator[](size_type idx) const
    {
        OCTOPUS_ASSERT_FMT_MSG(idx < data_.size(),
            "index (%1%) is larger than the data size (%2%)",
            idx % data_.size()); 
        return data_[idx];
    }
};

template <typename T0, typename T1>
bool same_dimensions(
    vector3d<T0> const& a
  , vector3d<T1> const& b
    )
{
    return (a.size() == b.size()) // Necessary but not sufficient.
        && (a.x_length() == b.x_length())
        && (a.y_length() == b.y_length())
        && (a.z_length() == b.z_length())
        ;
}

#define OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(OP)                     \
    template <typename T0, typename T1>                                       \
    vector3d<T0>& operator OP(                                                \
        vector3d<T0>& a                                                       \
      , T1 b                                                                  \
        )                                                                     \
    {                                                                         \
        using namespace octopus::operators;                                   \
        a.data_ OP b;                                                         \
        return a;                                                             \
    }                                                                         \
                                                                              \
    template <typename T>                                                     \
    vector3d<T>& operator OP(                                                 \
        vector3d<T>& a                                                        \
      , vector3d<T> const& b                                                  \
        )                                                                     \
    {                                                                         \
        using namespace octopus::operators;                                   \
        OCTOPUS_ASSERT(same_dimensions(a, b));                                \
        a.data_ OP b;                                                         \
        return a;                                                             \
    }                                                                         \
    /**/

OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(-=)
OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(+=)
OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(*=)
OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(/=)

#undef OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR

#define OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(OP)              \
    template <typename T0, typename T1>                     \
    vector3d<T0> operator OP(                               \
        vector3d<T0> const& a                               \
      , T1 b                                                \
        )                                                   \
    {                                                       \
        vector3d<T0> tmp(a);                                \
        return (tmp BOOST_PP_CAT(OP, =) b);                 \
    }                                                       \
                                                            \
    template <typename T>                                   \
    vector3d<T> operator OP(                                \
        vector3d<T> const& a                                \
      , vector3d<T> const& b                                \
        )                                                   \
    {                                                       \
        vector3d<T> tmp(a);                                 \
        return (tmp BOOST_PP_CAT(OP, =) b);                 \
    }                                                       \
    /**/

OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(-)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(+)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(*)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(/)

#undef OCTOPUS_DEFINE_ARITHMETIC_OPERATOR

}

#endif // OCTOPUS_820FC460_190E_4E57_8CB7_CB0880FC3E28

