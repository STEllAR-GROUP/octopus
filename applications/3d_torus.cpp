////////////////////////////////////////////////////////////////////////////////
//  Coypright (c) 2012 Zach Byerly 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/science.hpp>

struct enforce_outflow : octopus::trivial_serialization
{
    typedef void result_type;

    result_type operator()(boost::array<double, 3> x, octopus::face f) const
    {
        // IMPLEMENT
    } 
};

void octopus_define_problem(octopus::science_table& sci)
{
    ///////////////////////////////////////////////////////////////////////////
    sci.state_size = 1;

    sci.rho = octopus::get_value_from_state<0>();

    ///////////////////////////////////////////////////////////////////////////
    sci.enforce_outflow = enforce_outflow();
}


