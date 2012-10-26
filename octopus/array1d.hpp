////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_2B645EED_85D3_4967_89C9_0905C684BE8C)
#define OCTOPUS_2B645EED_85D3_4967_89C9_0905C684BE8C

#include <octopus/assert.hpp>

#include <boost/cstdint.hpp>
#include <boost/move/move.hpp>
#include <boost/serialization/access.hpp>

namespace octopus
{

template <typename T, boost::uint64_t Size>
struct array1d
{
    typedef boost::uint64_t size_type;

  private:
    T data_ [Size];

    BOOST_COPYABLE_AND_MOVABLE(array1d);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & boost::serialization::make_array(data_, Size);     
    }

  public:
    array1d()
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = boost::move(T());
    }

    array1d(array1d const& other)
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = other.data_[i];
    }

    array1d(BOOST_RV_REF(array1d) other)
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = boost::move(other.data_[i]);
    }

    array1d& operator=(BOOST_COPY_ASSIGN_REF(array1d) other)
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = other.data_[i];
        return *this;
    }

    array1d& operator=(BOOST_RV_REF(array1d) other)
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = boost::move(other.data_[i]);
        return *this;
    }

    array1d& operator=(T const& value)
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = value;
        return *this;
    }

    size_type size() const
    {
        return Size;
    } 

    T& operator()(size_type c)
    {
        OCTOPUS_ASSERT_FMT_MSG(c < Size,
            "coordinate (%1%) is larger than the length (%2%)",
            c % Size);  
        return data_[c];
    }

    T const& operator()(size_type c) const
    {
        OCTOPUS_ASSERT_FMT_MSG(c < Size,
            "coordinate (%1%) is larger than the length (%2%)",
            c % Size);  
        return data_[c];
    }

    T& operator[](size_type c)
    {
        return (*this)(c);
    }

    T const& operator[](size_type c) const
    {
        return (*this)(c);
    }
};

#define OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(OP)                     \
    template <typename T0, boost::uint64_t Size, typename T1>                 \
    array1d<T0, Size>& operator OP(                                           \
        array1d<T0, Size>& a                                                  \
      , T1 b                                                                  \
        )                                                                     \
    {                                                                         \
        a.data_ OP b;                                                         \
        return a;                                                             \
    }                                                                         \
                                                                              \
    template <typename T, boost::uint64_t Size>                               \
    array1d<T, Size>& operator OP(                                            \
        array1d<T, Size>& a                                                   \
      , array1d<T, Size> const& b                                             \
        )                                                                     \
    {                                                                         \
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

#define OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(OP)                  \
    template <typename T0, boost::uint64_t Size, typename T1>   \
    array1d<T0, Size> operator OP(                              \
        array1d<T0, Size> const& a                              \
      , T1 b                                                    \
        )                                                       \
    {                                                           \
        array1d<T0, Size> tmp(a);                               \
        return (tmp BOOST_PP_CAT(OP, =) b);                     \
    }                                                           \
                                                                \
    template <typename T0, boost::uint64_t Size, typename T1>   \
    array1d<T0, Size> operator OP(                              \
        array1d<T0, Size>&& a                                   \
      , T1 b                                                    \
        )                                                       \
    {                                                           \
        array1d<T0, Size> tmp(a);                               \
        return (tmp BOOST_PP_CAT(OP, =) b);                     \
    }                                                           \
                                                                \
    template <typename T, boost::uint64_t Size>                 \
    array1d<T, Size> operator OP(                               \
        array1d<T, Size> const& a                               \
      , array1d<T, Size> const& b                               \
        )                                                       \
    {                                                           \
        array1d<T, Size> tmp(a);                                \
        return (tmp BOOST_PP_CAT(OP, =) b);                     \
    }                                                           \
                                                                \
    template <typename T, boost::uint64_t Size>                 \
    array1d<T, Size> operator OP(                               \
        array1d<T, Size>&& a                                    \
      , array1d<T, Size> const& b                               \
        )                                                       \
    {                                                           \
        array1d<T, Size> tmp(a);                                \
        return (tmp BOOST_PP_CAT(OP, =) b);                     \
    }                                                           \
    /**/

OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(-)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(+)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(*)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(/)

#undef OCTOPUS_DEFINE_ARITHMETIC_OPERATOR

}

#endif // OCTOPUS_2B645EED_85D3_4967_89C9_0905C684BE8C

