////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Dominic Marcello
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_0091D3AA_2B0C_4069_80E5_9EA08E57C45F)
#define OCTOPUS_0091D3AA_2B0C_4069_80E5_9EA08E57C45F

#include <octopus/config.hpp>
#include <octopus/state.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT minmod_reconstruction
{
  private:
    double theta_;

  public:
    minmod_reconstruction(double theta = 1.3) : theta_(theta) {}

    enum { ghost_zone_width = 2 };

    void operator()(
        std::vector<state> const& q0
      , std::vector<state>& ql
      , std::vector<state>& qr
        ) const;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & theta_;
    }
};

}

#endif // OCTOPUS_0091D3AA_2B0C_4069_80E5_9EA08E57C45F

