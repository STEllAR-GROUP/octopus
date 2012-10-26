////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_F5DDF845_26AB_417F_9AC2_2824296CC18A)
#define OCTOPUS_F5DDF845_26AB_417F_9AC2_2824296CC18A

#include <octopus/assert.hpp>

#include <boost/preprocessor/cat.hpp>
#include <boost/move/move.hpp>

#include <vector>

namespace octopus { namespace operators
{

#define OCTOPUS_DEFINE_ARITHMETIC_ASSIGNMENT_OPERATOR(OP)                     \
    template <typename T0, typename T1>                                       \
    std::vector<T0>& operator OP(                                             \
        std::vector<T0>& a                                                    \
      , T1 b                                                                  \
        )                                                                     \
    {                                                                         \
        for (std::size_t i = 0; i < a.size(); ++i)                            \
            a[i] OP b;                                                        \
        return a;                                                             \
    }                                                                         \
                                                                              \
    template <typename T>                                                     \
    std::vector<T>& operator OP(                                              \
        std::vector<T>& a                                                     \
      , std::vector<T> const& b                                               \
        )                                                                     \
    {                                                                         \
        OCTOPUS_ASSERT_MSG(a.size() == b.size(),                              \
            "dimensions do not match");                                       \
        for (std::size_t i = 0; i < a.size(); ++i)                            \
            a[i] OP b[i];                                                     \
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
    std::vector<T0> operator OP(                            \
        std::vector<T0> const& a                            \
      , T1 b                                                \
        )                                                   \
    {                                                       \
        std::vector<T0> tmp(a);                             \
        return (tmp BOOST_PP_CAT(OP, =) b);                 \
    }                                                       \
                                                            \
    template <typename T0, typename T1>                     \
    std::vector<T0> operator OP(                            \
        BOOST_RV_REF(std::vector<T0>) a                     \
      , T1 b                                                \
        )                                                   \
    {                                                       \
        std::vector<T0> tmp(a);                             \
        return (tmp BOOST_PP_CAT(OP, =) b);                 \
    }                                                       \
                                                            \
    template <typename T>                                   \
    std::vector<T> operator OP(                             \
        std::vector<T> const& a                             \
      , std::vector<T> const& b                             \
        )                                                   \
    {                                                       \
        std::vector<T> tmp(a);                              \
        return (tmp BOOST_PP_CAT(OP, =) b);                 \
    }                                                       \
                                                            \
    template <typename T>                                   \
    std::vector<T> operator OP(                             \
        BOOST_RV_REF(std::vector<T>) a                      \
      , std::vector<T> const& b                             \
        )                                                   \
    {                                                       \
        std::vector<T> tmp(a);                              \
        return (tmp BOOST_PP_CAT(OP, =) b);                 \
    }                                                       \
    /**/

OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(-)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(+)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(*)
OCTOPUS_DEFINE_ARITHMETIC_OPERATOR(/)

#undef OCTOPUS_DEFINE_ARITHMETIC_OPERATOR

template <typename T>
std::vector<T> operator-(
    std::vector<T> const& a
    ) 
{
    std::vector<T> tmp; 
    for (std::size_t i = 0; i < a.size(); ++i)
        tmp[i] = T(0) - a[i];
    return tmp; 
}

template <typename T>
std::vector<T> operator-(
    BOOST_RV_REF(std::vector<T>) a
    ) 
{
    std::vector<T> tmp(a); 
    for (std::size_t i = 0; i < a.size(); ++i)
        tmp[i] = T(0) - tmp[i];
    return tmp; 
}

template <typename T>
std::vector<T> operator+(
    std::vector<T> const& a
    ) 
{
    return a;
}

template <typename T>
std::vector<T> operator+(
    BOOST_RV_REF(std::vector<T>) a
    ) 
{
    return a;
}


}}

#endif // OCTOPUS_F5DDF845_26AB_417F_9AC2_2824296CC18A

