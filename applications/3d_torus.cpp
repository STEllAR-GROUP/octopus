////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

// http://www.vistrails.org/index.php/User:Tohline/Apps/PapaloizouPringleTori

#include <octopus/driver.hpp>
#include <octopus/science.hpp>
#include <octopus/octree/checkout_state.hpp>

/// Mass density
inline double&       rho(std::vector<double>& u)       { return u.at(0); }
inline double const& rho(std::vector<double> const& u) { return u.at(0); }

/// Momentum density (X-axis)
inline double&       sx(std::vector<double>& u)       { return u.at(1); }
inline double const& sx(std::vector<double> const& u) { return u.at(1); }

/// Momentum density (Y-axis)
inline double&       sy(std::vector<double>& u)       { return u.at(2); }
inline double const& sy(std::vector<double> const& u) { return u.at(2); }

/// Momentum density (Z-axis)
inline double&       sz(std::vector<double>& u)       { return u.at(3); }
inline double const& sz(std::vector<double> const& u) { return u.at(3); }

/// Total energy (Zach was not positive) 
inline double&       et(std::vector<double>& u)       { return u.at(4); }
inline double const& et(std::vector<double> const& u) { return u.at(4); }

/// Entropy tracer
inline double&       tau(std::vector<double>& u)       { return u.at(5); }
inline double const& tau(std::vector<double> const& u) { return u.at(5); }

// Initialization kernel.
struct initialize : octopus::trivial_serialization
{
    void operator()(octopus::octree_server& e) const
    {
        octopus::checkout_state U(e, octopus::checkout_for_init);

        using std::pow;
        using std::sqrt;

        // Polytropic index.
        double const gamma = 2.0; // EULER_GAMMA

        // Polytropic constant.
        double const kappa = 1.0;

        double const ei0 = 1.0;
        double const tau0 = pow(ei0, 1.0 / gamma);
        double const rho1 = 1.0e-10;
        double const ei1 = 1.0e-10;
        double const tau1 = pow(ei1, 1.0 / gamma);
    
        double const G = 1.0;

        // Mass of the central object.
        double const M_c = 2e-2; // M_C
    
        double const eps = 0.4;
        double const R_outer = 1.0747e-4;
        double const R_inner = R_outer*(1.0 - eps)/(1.0 + eps);

        double const h = sqrt(2.0*G*M_c*R_inner*R_outer/(R_inner + R_outer));

        double const C = 0.5*pow(h/R_inner, 2) - G*M_c/R_inner;    
   
 
        for (boost::uint64_t i = 0; i < U.x_length(); ++i)
        {
            for (boost::uint64_t j = 0; j < U.y_length(); ++j)
            {
                for (boost::uint64_t k = 0; k < U.z_length(); ++k)
                {
                    double const x = e.xc(i);
                    double const y = e.yc(j);
                    double const z = e.zc(k);
                  
                    // Cylindrical R.  
                    double const r = sqrt(pow(x, 2) + pow(y, 2));
    
                    if ((R_inner <= r) && (R_outer >= r))
                    {
                        double const z_max =
                            sqrt(pow(G*M_c/(0.5*pow(h/r, 2) - C), 2) - r*r);

                        if (z <= z_max)
                        {
                            double const rho_here =
                                  (0.5/kappa)
                                * (C + G*M_c/sqrt(r*r + z*z) - 0.5*pow(h/r, 2));

                            rho(U(i, j, k)) = rho_here;
                            sx(U(i, j, k))  = -y*rho_here*h/pow(r, 2);
                            sy(U(i, j, k))  = x*rho_here*h/pow(r, 2);
                            et(U(i, j, k))  = ei0;
                            tau(U(i, j, k)) = tau0;
                        }
                    }
                    
                    else
                    {
                        rho(U(i, j, k)) = rho1;
                        sx(U(i, j, k))  = 0.0; 
                        sy(U(i, j, k))  = 0.0;
                        et(U(i, j, k))  = ei1;  
                        tau(U(i, j, k)) = tau1;
                    }

                    sz(U(i, j, k)) = 0.0;
                }
            }
        }
    }
};

struct enforce_outflow : octopus::trivial_serialization
{
    typedef void result_type;

    result_type operator()(octopus::face f, boost::array<double, 3> x) const
    {
        // IMPLEMENT
    } 
};

void octopus_define_problem(octopus::science_table& sci)
{
    sci.state_size = 6;

    sci.enforce_outflow = enforce_outflow();

    sci.initialize = initialize();
}


