////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_C509AE65_1FE4_4FAD_96B5_BB75605EB37C)
#define OCTOPUS_C509AE65_1FE4_4FAD_96B5_BB75605EB37C

#include <hpx/exception.hpp>

#include <octopus/config.hpp>

#include <boost/current_function.hpp>

#if OCTOPUS_VERIFY
    #define OCTOPUS_ASSERT(expr) ((expr)                                \
      ? ((void)0)                                                       \
      : ::boost::assertion_failed                                       \
            (#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))

    #define OCTOPUS_ASSERT_MSG(expr, msg) ((expr)                       \
      ? ((void)0)                                                       \
      : ::boost::assertion_failed_msg                                   \
            (#expr, msg, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))
#else
    #define OCTOPUS_ASSERT(expr)
    #define OCTOPUS_ASSERT_MSG(expr)
#endif

#endif // OCTOPUS_C509AE65_1FE4_4FAD_96B5_BB75605EB37C

