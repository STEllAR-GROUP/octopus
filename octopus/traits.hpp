////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_B5F2FEC8_0541_4FBA_BCD3_268E2C0DA912)
#define OCTOPUS_B5F2FEC8_0541_4FBA_BCD3_268E2C0DA912

#include <hpx/traits.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_pod.hpp>

namespace octopus 
{

template <typename T, typename Enable = void>
struct parameter_type
{
    typedef T const& type;
};

template <typename T>
struct parameter_type<T,
    typename boost::enable_if<
        boost::mpl::or_<
            boost::is_arithmetic<T>
          , boost::is_pod<T>
          , boost::serialization::is_bitwise_serializable<T>
        >
      , void
    >::type
> {
    typedef T type;
};

template <typename T, typename Enable = void>
struct is_bool : boost::mpl::false_ {};

template <>
struct is_bool<bool> : boost::mpl::true_ {};

template <typename T>
struct proxied_type
{
    typedef T type;
};

}

#endif // OCTOPUS_B5F2FEC8_0541_4FBA_BCD3_268E2C0DA912

