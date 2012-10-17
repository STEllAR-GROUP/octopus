////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_747228A8_81BD_4AEB_A537_842E687DA4EE)
#define OCTOPUS_747228A8_81BD_4AEB_A537_842E687DA4EE

#include <octopus/assert.hpp>

#include <boost/preprocessor/cat.hpp>
#include <boost/array.hpp>

namespace octopus { namespace operators
{

#define OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(OP)                     \
    template <typename T0, std::size_t Size, typename T1>                     \
    boost::array<T0, Size>& operator OP(                                      \
        boost::array<T0, Size>& a                                             \
      , T1 b                                                                  \
        )                                                                     \
    {                                                                         \
        for (std::size_t i = 0; i < Size; ++i)                                \
            a[i] OP b;                                                        \
        return a;                                                             \
    }                                                                         \
                                                                              \
    template <typename T, std::size_t Size>                                   \
    boost::array<T, Size>& operator OP(                                       \
        boost::array<T, Size>& a                                              \
      , boost::array<T, Size> const& b                                        \
        )                                                                     \
    {                                                                         \
        for (std::size_t i = 0; i < Size; ++i)                                \
            a[i] OP b[i];                                                     \
        return a;                                                             \
    }                                                                         \
    /**/

OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(-=)
OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(+=)
OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(*=)
OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(/=)

#undef OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR

#define OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(OP)                     \
    template <typename T0, std::size_t Size, typename T1>          \
    boost::array<T0, Size> operator OP(                            \
        boost::array<T0, Size> const& a                            \
      , T1 b                                                       \
        )                                                          \
    {                                                              \
        boost::array<T0, Size> tmp = a;                            \
        return (tmp BOOST_PP_CAT(OP, =) b);                        \
    }                                                              \
                                                                   \
    template <typename T, std::size_t Size>                        \
    boost::array<T, Size> operator OP(                             \
        boost::array<T, Size> const& a                             \
      , boost::array<T, Size> const& b                             \
        )                                                          \
    {                                                              \
        boost::array<T, Size> tmp = a;                             \
        return (tmp BOOST_PP_CAT(OP, =) b);                        \
    }                                                              \
    /**/

OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(-)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(+)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(*)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(/)

#undef OCTOPUS_DEFINE_ARITHMETIC_OPERATOR

template <typename T, std::size_t Size>
boost::array<T, Size> operator-(
    boost::array<T, Size> const& a
    ) 
{
    boost::array<T, Size> tmp; 
    for (std::size_t i = 0; i < Size; ++i)
        tmp[i] = T(0) - a[i];
    return tmp; 
}

template <typename T, std::size_t Size>
boost::array<T, Size> operator+(
    boost::array<T, Size> const& a
    ) 
{
    return a;
}

template <typename T0, std::size_t Size, typename T1>
boost::array<T0, Size>& operator<<=(
    boost::array<T0, Size>& a
  , T1 b
    )
{
    for (std::size_t i = 0; i < Size; ++i)
        a[i] <<= b;
    return a;
}

template <typename T0, std::size_t Size, typename T1>
boost::array<T0, Size>& operator>>=(
    boost::array<T0, Size>& a
  , T1 b
    )
{
    for (std::size_t i = 0; i < Size; ++i)
        a[i] >>= b;
    return a;
}

}}

#endif // OCTOPUS_747228A8_81BD_4AEB_A537_842E687DA4EE

