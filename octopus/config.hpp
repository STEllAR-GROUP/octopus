////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2007-2012 Hartmut Kaiser
//  Copyright (c) 2011-2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_974D87CF_23F9_41CC_9705_2A4CDB3339DC)
#define OCTOPUS_974D87CF_23F9_41CC_9705_2A4CDB3339DC

#include <hpx/config.hpp>

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

#if defined(__GNUC__)
    #define OCTOPUS_WEAK __attribute__((weak))
#else
    #define OCTOPUS_WEAK
#endif

#if defined(_DEBUG) && !defined(DEBUG)
    #define DEBUG
#endif

#if defined(DEBUG) && !defined(OCTOPUS_DEBUG)
    #define OCTOPUS_DEBUG 1
#endif

#if !defined(OCTOPUS_ENABLE_VERIFICATION)
    #define OCTOPUS_ENABLE_VERIFICATION 0
#endif

#if !defined(OCTOPUS_STATE_SIZE)
    #define OCTOPUS_STATE_SIZE 9
#endif

#if BOOST_VERSION < 105300
    #if defined(BOOST_NO_RVALUE_REFERENCES)
        #define BOOST_COPY_ASSIGN_REF_3_TEMPL_ARGS \
                BOOST_MOVE_COPY_ASSIGN_REF_3_TEMPL_ARGS
    #else
        #define BOOST_COPY_ASSIGN_REF_3_TEMPL_ARGS(TYPE, ARG1, ARG2, ARG3) \
            const TYPE<ARG1, ARG2, ARG3>&
    #endif
#endif

#endif // OCTOPUS_974D87CF_23F9_41CC_9705_2A4CDB3339DC

