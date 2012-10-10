////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/util/lightweight_test.hpp>

#include <octopus/driver.hpp>
#include <octopus/vector3d.hpp>

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::vector3d<boost::uint64_t> v;

    return hpx::util::report_errors();
}

