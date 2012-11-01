////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_286AA42C_449B_4478_B903_EE7F7177FFE8)
#define OCTOPUS_286AA42C_449B_4478_B903_EE7F7177FFE8

#include <boost/filesystem.hpp>

namespace octopus
{

// How is this (or operator/) not a thing in Boost.Filesystem?
inline std::string join_paths(
    std::string const& a
  , std::string const& b
    )
{
    namespace fs = boost::filesystem;

    fs::path aa(a);
    aa /= b;

    return aa.string();
}

inline std::string join_paths(
    std::string const& a
  , std::string const& b
  , std::string const& c
    )
{
    namespace fs = boost::filesystem;

    fs::path aa(a);
    aa /= b;

    return aa.string();
}

inline std::string current_path()
{
    return boost::filesystem::current_path().string();
}

}

#endif // OCTOPUS_286AA42C_449B_4478_B903_EE7F7177FFE8

