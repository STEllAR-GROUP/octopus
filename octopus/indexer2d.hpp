////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_0FC14503_A70D_41C0_A590_60F982BAFCDA)
#define OCTOPUS_0FC14503_A70D_41C0_A590_60F982BAFCDA

#include <octopus/assert.hpp>

namespace octopus
{

/// Utility from Dominic's original code, I have retained this to ease porting
/// and keep the two codes as similar as possible. Step = 1 is the equivalent of
/// the Indexer2d class from the original code; Step = 2 is the equivalent of
/// Indexer2d_by2. 
template <
    std::size_t Step = 1
>
struct indexer2d
{
  private:
    std::size_t const xstart_;
    std::size_t const ystart_;
    std::size_t const xspan_;
    std::size_t const yspan_;

  public:
    std::size_t const maximum;

    indexer2d(
        std::size_t xstart,
        std::size_t xstop,
        std::size_t ystart,
        std::size_t ystop
        )
      : xstart_(xstart)
      , ystart_(ystart)
      , xspan_(xstop - xstart + 1)
      , yspan_(ystop - ystart + 1)
      , maximum((ystop - ystart + 1) * (xstop - xstart + 1) / (Step * Step) - 1)
        // NOTE: Not sure about the division by (Step * Step). In the original
        // Indexer2d (Step = 1), there is no division there. In Indexer2d_by2
        // (Step = 2), the denominator is 4. I'm fairly certain that
        // (Step * Step) works for Step > 2, but this is not tested.
    {
        OCTOPUS_ASSERT(0 < xspan_);
        OCTOPUS_ASSERT(0 < yspan_);
    }

    std::size_t x(std::size_t idx) const
    {
        OCTOPUS_ASSERT(idx >= 0);
        OCTOPUS_ASSERT(idx <= maximum);
        return ((Step * idx) % xspan_) + xstart_; 
    }

    std::size_t y(std::size_t idx) const
    {
        OCTOPUS_ASSERT(idx >= 0);
        OCTOPUS_ASSERT(idx <= maximum);
        return Step * ((Step * idx) / xspan_) + ystart_;
    }
};

}

#endif // OCTOPUS_0FC14503_A70D_41C0_A590_60F982BAFCDA

