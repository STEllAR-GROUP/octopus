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

#if OCTOPUS_ENABLE_VERIFICATION
    #include <boost/current_function.hpp>
    #include <boost/format.hpp>

    #define OCTOPUS_ASSERT(expr) ((expr)                                \
      ? ((void)0)                                                       \
      : ::boost::assertion_failed                                       \
            (#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))

    #define OCTOPUS_ASSERT_MSG(expr, msg) ((expr)                       \
      ? ((void)0)                                                       \
      : ::boost::assertion_failed_msg                                   \
            (#expr, msg, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))

    #define OCTOPUS_ASSERT_FMT_MSG(expr, fmt, args) ((expr)             \
      ? ((void)0)                                                       \
      : ::boost::assertion_failed_msg                                   \
            (#expr, boost::str(boost::format(fmt) % args).c_str(),      \
                BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))
#else
    #define OCTOPUS_ASSERT(expr)
    #define OCTOPUS_ASSERT_MSG(expr)
    #define OCTOPUS_ASSERT_FMT_MSG(expr, fmt, args)
#endif

#if OCTOPUS_ENABLE_TEST_IN_PLACE
    #include <boost/current_function.hpp>

    #define OCTOPUS_TEST_IN_PLACE(expr) ((expr)                         \
      ? ((void)0)                                                       \
      : ::boost::assertion_failed                                       \
            (#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__))
#else
    #define OCTOPUS_TIP(expr)
#endif

#endif // OCTOPUS_C509AE65_1FE4_4FAD_96B5_BB75605EB37C

