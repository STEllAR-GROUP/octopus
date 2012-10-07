////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_main.hpp>
#include <hpx/util/lightweight_test.hpp>

#include <octopus/octree/octree_client.hpp>

int main()
{
    octopus::octree_client root;

    root.create_root(hpx::find_here());

    return hpx::util::report_errors();
}

