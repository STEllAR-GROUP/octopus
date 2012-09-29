////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICESizeSE_1_0.txt or copy at http://www.boost.org/LICESizeSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_4765E5BC_E984_4DF3_B1EE_84B1B8D3C4A7)
#define OCTOPUS_4765E5BC_E984_4DF3_B1EE_84B1B8D3C4A7

#include <boost/array.hpp>
#include <boost/serialization/array.hpp>

#include <cmath>

namespace octopus
{

template <typename T, std::size_t Size>
struct array1d : boost::array<T, Size>
{
    typedef boost::array<T, Size> base_type;

    array1d()
    {
        *this = T();
    }

    array1d(array1d const& v)
    {
        *this = v;
    }

    array1d(T const& v)
    {
        *this = v;
    }

    array1d& operator=(array1d const& v)
    {
        base_type::operator=(v);
        return *this;
    }

    array1d& operator=(T const& v)
    {
        for (std::size_t i = 0; i < Size; ++i) 
            (*this)[i] = v;
        return *this;
    }

    T dot(array1d const& v) const
    {
        T d = 0;
        for (std::size_t i = 0; i < Size; ++i) 
            d += (*this)[i] * v[i];
        return d;
    }

    T mag() const
    {
        return std::sqrt(this->dot(*this));
    }

    array1d& operator/=(T const& t)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] /= t;
        return *this;
    }

    array1d& operator*=(T const& t)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] *= t;
        return *this;
    }

    array1d& operator-=(T const& t)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] -= t;
        return *this;
    }

    array1d& operator+=(T const& t)
    {
        for (std::size_t i = 0; i < Size; ++i) 
            (*this)[i] += t;
        return *this;
    }

    array1d& operator+=(array1d const& v)
    {
        for (std::size_t i = 0; i < Size; ++i) 
            (*this)[i] += v[i];
        return *this;
    }

    array1d& operator-=(array1d const& v)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] -= v[i];
        return *this;
    }

    array1d& operator*=(array1d const& v)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] *= v[i];
        return *this;
    }

    array1d& operator/=(array1d const& v)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] /= v[i];
        return *this;
    }

    array1d operator+(array1d const& v) const
    {
        return (array1d(*this) += v);
    }

    array1d operator-(array1d const& v) const
    {
        return (array1d(*this) -= v);
    }

    array1d operator*(array1d const& v) const
    {
        return (array1d(*this) *= v);
    }

    array1d operator/(array1d const& v) const
    {
        return (array1d(*this) /= v);
    }

    array1d operator+(T const& v) const
    {
        return (array1d(*this) += v);
    }

    array1d operator-(T const& v) const
    {
        return (array1d(*this) -= v);
    }

    array1d operator*(T const& v) const
    {
        return (array1d(*this) *= v);
    }

    array1d operator/(T const& v) const
    {
        return (array1d(*this) /= v);
    }

    array1d operator-() const
    {
        return (array1d(T(0)) -= *this);
    }

    array1d operator+() const
    {
        return array1d(*this);
    }

    array1d& operator<<=(T const& v)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] <<= v;
        return *this;
    }

    array1d& operator>>=(T const& v)
    {
        for (std::size_t i = 0; i < Size; ++i)
            (*this)[i] >>= v;
        return *this;
    }
};

}

#endif // OCTOPUS_4765E5BC_E984_4DF3_B1EE_84B1B8D3C4A7

