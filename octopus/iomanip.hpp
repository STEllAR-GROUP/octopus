////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_212F0E1D_32DD_43AD_A3AE_E14565F2AB2D)
#define OCTOPUS_212F0E1D_32DD_43AD_A3AE_E14565F2AB2D

#include <ostream>

namespace octopus
{

template <typename T>
inline T highbit()
{
    return (((T)(-1)) >> 1) + 1;
}

template <typename T>
struct binary_formatter
{
  private:
    T value_;

  public: 
    binary_formatter(T value) : value_(value) {}

    binary_formatter(binary_formatter const& other) : value_(other.value_) {}

    friend std::ostream& operator<<(std::ostream& os, binary_formatter f)
    {
        for (T bit = highbit<T>(); bit; bit >>= 1)
            os << ((f.value_ & bit) ? '1' : '0');
        return os;
    }
};

template <typename T>
inline binary_formatter<T> binary(T value)
{
    return binary_formatter<T>(value);
}

}

#endif // OCTOPUS_212F0E1D_32DD_43AD_A3AE_E14565F2AB2D

