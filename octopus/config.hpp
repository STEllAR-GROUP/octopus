////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2007-2012 Hartmut Kaiser
//  Copyright (c) 2011-2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_974D87CF_23F9_41CC_9705_2A4CDB3339DC)
#define OCTOPUS_974D87CF_23F9_41CC_9705_2A4CDB3339DC

#include <boost/config.hpp>

#if defined(BOOST_MSVC) 
    #define OCTOPUS_SYMBOL_EXPORT __declspec(dllexport)
    #define OCTOPUS_SYMBOL_IMPORT __declspec(dllimport)
#else
    #define OCTOPUS_SYMBOL_EXPORT __attribute__((visibility("default")))
    #define OCTOPUS_SYMBOL_IMPORT __attribute__((visibility("default")))
#endif

#if defined(OCTOPUS_EXPORTS)
    #define OCTOPUS_EXPORT OCTOPUS_SYMBOL_EXPORT
#else
    #define OCTOPUS_EXPORT OCTOPUS_SYMBOL_IMPORT
#endif

#if defined(OCTOPUS_EXPORTS)
    #define OCTOPUS_EXPORT OCTOPUS_SYMBOL_EXPORT
#else
    #define OCTOPUS_EXPORT OCTOPUS_SYMBOL_IMPORT
#endif

#if defined(_DEBUG) && !defined(DEBUG)
    #define DEBUG
#endif

#if defined(DEBUG) && !defined(OCTOPUS_DEBUG)
    #define OCTOPUS_DEBUG 1
#endif

#if !defined(OCTOPUS_ENABLE_VERIFICATION)
    #define OCTOPUS_ENABLE_VERIFICATION 1
#endif

#if !defined(OCTOPUS_ENABLE_TEST_IN_PLACE)
    #define OCTOPUS_ENABLE_TEST_IN_PLACE 1
#endif

#endif // OCTOPUS_974D87CF_23F9_41CC_9705_2A4CDB3339DC

