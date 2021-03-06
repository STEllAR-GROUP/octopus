////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_754BFB65_F28B_4F97_83FF_69C5DEAB500F)
#define OCTOPUS_754BFB65_F28B_4F97_83FF_69C5DEAB500F

#include <octopus/assert.hpp>
#include <octopus/trivial_serialization.hpp>
#include <octopus/state.hpp>

#include <boost/cstdint.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace octopus
{

template <typename T>
inline bool compare_real(T x, T y)
{
    T const epsilon = std::numeric_limits<T>::epsilon();
 
    if ((x + epsilon >= y) && (x - epsilon <= y))
        return true;
    else
        return false;
}

template <typename T>
inline bool compare_real(T x, T y, T epsilon)
{
    if ((x + epsilon >= y) && (x - epsilon <= y))
        return true;
    else
        return false;
}

template <typename T>
inline boost::uint64_t hex_real(T x)
{
    boost::uint64_t* hex = (boost::uint64_t*) &x;
    return *hex;
}

inline double maximum(double a, double b)
{
    return (std::max)(a, b);
}

inline double maximum(double a, double b, double c)
{
    return (std::max)(a, (std::max)(b, c));
}

inline double minimum(double a, double b)
{
    return (std::min)(a, b);
}

inline double minimum(double a, double b, double c)
{
    return (std::min)(a, (std::min)(b, c));
}

// Serializable minimum
struct minimum_functor : octopus::trivial_serialization
{
    template <typename T>
    T const& operator()(T const& a, T const& b) const
    {
        return (std::min)(a, b);
    }
};

// Serializable maximum
struct maximum_functor : octopus::trivial_serialization
{
    template <typename T>
    T const& operator()(T const& a, T const& b) const
    {
        return (std::max)(a, b);
    }
};

template <typename T>
inline T sign(T a)
{
    if (a > 0) 
        return 1;
    else if (a < 0) 
        return -1;
    else 
        return 0;
}

///////////////////////////////////////////////////////////////////////////////
inline double minmod(double a, double b)
{
    return 0.5 * (sign(a) + sign(b)) * (std::min)(std::fabs(a), std::fabs(b));
}

inline double minmod(double a, double b, double c)
{
    return minmod(a, minmod(b, c));
}

inline double minmod_theta(double a, double b, double theta)
{
    return minmod(theta * a, theta * b, 0.5 * (a + b));
}

///////////////////////////////////////////////////////////////////////////////
template <std::size_t Size, typename Rep0, typename Rep1>
inline array<double, Size> minmod(
    array<double, Size, Rep0> const& v1
  , array<double, Size, Rep1> const& v2
    )
{
    //OCTOPUS_ASSERT(v1.size() == v2.size());

    array<double, Size> mm;
    //mm.reserve(v1.size());

    for (std::size_t i = 0; i < v1.size(); ++i) 
        mm[i] = minmod(v1[i], v2[i]);

    return mm;
}

template <std::size_t Size, typename Rep0, typename Rep1, typename Rep2>
inline array<double, Size> minmod(
    array<double, Size, Rep0> const& v1
  , array<double, Size, Rep1> const& v2
  , array<double, Size, Rep2> const& v3)
{
    //OCTOPUS_ASSERT(v1.size() == v2.size());
    //OCTOPUS_ASSERT(v1.size() == v3.size());

    array<double, Size> mm;
    //mm.reserve(v1.size());

    for (std::size_t i = 0; i < v1.size(); ++i) 
        mm[i] = minmod(v1[i], v2[i], v3[i]);

    return mm;
}

template <std::size_t Size, typename Rep0, typename Rep1>
inline array<double, Size> minmod_theta(
    array<double, Size, Rep0> const& v1
  , array<double, Size, Rep1> const& v2
  , double theta
    )
{
    OCTOPUS_ASSERT(v1.size() == v2.size());

    array<double, Size> mm;
    //mm.reserve(v1.size());

    for (std::size_t i = 0; i < v1.size(); ++i) 
        mm[i] = minmod(
            theta * v1[i]
          , theta * v2[i]
          , 0.5 * (v1[i] + v2[i]));

    return mm;
}

}

#endif // OCTOPUS_754BFB65_F28B_4F97_83FF_69C5DEAB500F

