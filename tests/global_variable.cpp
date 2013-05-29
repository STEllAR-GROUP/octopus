////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2013 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <hpx/hpx_main.hpp>
#include <hpx/util/lightweight_test.hpp>

#include <octopus/global_variable.hpp>

OCTOPUS_GLOBAL_VARIABLE((int), test0);
OCTOPUS_GLOBAL_VARIABLE((int), test1, (17));

///////////////////////////////////////////////////////////////////////////////
int check_value0() { return test0; }
HPX_PLAIN_ACTION(check_value0, check_value0_action);

int check_value1() { return test1; }
HPX_PLAIN_ACTION(check_value1, check_value1_action);

///////////////////////////////////////////////////////////////////////////////
int main()
{
    test0 = 42;

    HPX_TEST_EQ(test0.get(), 42);

    {
        std::vector<hpx::id_type> targets = hpx::find_all_localities();

        std::vector<hpx::future<int> > futures;
        futures.reserve(targets.size()); 

        for (boost::uint64_t i = 0; i < targets.size(); ++i)
            futures.push_back(hpx::async<check_value0_action>(targets[i])); 

        hpx::lcos::wait(futures); 

        for (boost::uint64_t i = 0; i < targets.size(); ++i)
            HPX_TEST_EQ(futures[i].get(), 42);
    }

    HPX_TEST_EQ(test1.get(), 17);

    {
        std::vector<hpx::id_type> targets = hpx::find_all_localities();

        std::vector<hpx::future<int> > futures;
        futures.reserve(targets.size()); 

        for (boost::uint64_t i = 0; i < targets.size(); ++i)
            futures.push_back(hpx::async<check_value1_action>(targets[i])); 

        hpx::lcos::wait(futures); 

        for (boost::uint64_t i = 0; i < targets.size(); ++i)
            HPX_TEST_EQ(futures[i].get(), 17);
    }

    return hpx::util::report_errors();
}

