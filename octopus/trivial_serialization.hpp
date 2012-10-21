////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_F3666EE4_3AE3_4BFD_AF4B_9258F84DE803)
#define OCTOPUS_F3666EE4_3AE3_4BFD_AF4B_9258F84DE803

namespace octopus
{

/// Inherit from this class to implement a trivial serialize() method.
struct trivial_serialization
{
    // Trivial serialization support.
    template <typename Archive> void serialize(Archive&, unsigned int) {}
};

}

#endif // OCTOPUS_F3666EE4_3AE3_4BFD_AF4B_9258F84DE803

