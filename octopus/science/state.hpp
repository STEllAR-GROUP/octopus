////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_60B35C73_7FD4_47AB_8928_3AD5AC5FE5B5)
#define OCTOPUS_60B35C73_7FD4_47AB_8928_3AD5AC5FE5B5

#include <octopus/trivial_serialization.hpp>

namespace octopus
{

template <std::size_t Index>
struct get_value_from_state : trivial_serialization
{
    typedef double result_type;
 
    result_type operator()(
        std::vector<double> const& state 
        ) const
    {
        #if defined(OCTOPUS_VERIFY)
            return state.at(Index);
        #else
            return state[Index];
        #endif
    }
};

}

#endif // OCTOPUS_60B35C73_7FD4_47AB_8928_3AD5AC5FE5B5

