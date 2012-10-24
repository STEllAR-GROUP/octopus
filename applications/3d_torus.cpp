////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

// http://www.vistrails.org/index.php/User:Tohline/Apps/PapaloizouPringleTori

#include <octopus/driver.hpp>
#include <octopus/science.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/engine/ini.hpp>
#include <octopus/io/silo.hpp>

#include <octopus/operators/boost_array_arithmetic.hpp>

// FIXME: Proper configuration

// NOTE: Zach, explain what this is.    
double const G = 1.0; // BIG_G

// Mass of the central object.
double const M_C = 2e-2;

// Minimum value that rho is allowed to be.
double const rho_floor = 1.0e-20;

// Minimum value that internal energy is allowed to be. 
double const internal_energy_floor = 1.0e-20;  

// Polytropic index.
double const GAMMA = 2.0; // EULER_GAMMA

// Polytropic constant.
double const KAPPA = 1.0;

///////////////////////////////////////////////////////////////////////////////
/// Mass density
inline double&       rho(std::vector<double>& u)       { return u[0]; }
inline double const& rho(std::vector<double> const& u) { return u[0]; }

/// Momentum density (X-axis)
inline double&       momentum_x(std::vector<double>& u)       { return u[1]; }
inline double const& momentum_x(std::vector<double> const& u) { return u[1]; }

/// Momentum density (Y-axis)
inline double&       momentum_y(std::vector<double>& u)       { return u[2]; }
inline double const& momentum_y(std::vector<double> const& u) { return u[2]; }

/// Momentum density (Z-axis)
inline double&       momentum_z(std::vector<double>& u)       { return u[3]; }
inline double const& momentum_z(std::vector<double> const& u) { return u[3]; }

/// Total energy of the gas 
inline double&       total_energy(std::vector<double>& u)       { return u[4]; }
inline double const& total_energy(std::vector<double> const& u) { return u[4]; }

/// Entropy tracer
inline double&       tau(std::vector<double>& u)       { return u[5]; }
inline double const& tau(std::vector<double> const& u) { return u[5]; }

inline double kinetic_energy(std::vector<double> const& state)
{
    return 0.5 * ( momentum_x(state) * momentum_x(state)
                 + momentum_y(state) * momentum_y(state)
                 + momentum_z(state) * momentum_z(state)) / rho(state);
}

template <octopus::axis Axis>
inline double gravity(double x, double y, double z)
{
    using std::sqrt;

    double const r = sqrt(x*x + y*y + z*z);
    double const F = -G*M_C/(r*r);

    if (Axis == octopus::z_axis)
    {
        double const r_cyl = sqrt(x*x + y*y); 
        return F*(z/r_cyl);
    }

    return F*(y/r);
}

template <octopus::axis Axis>
inline double gravity(boost::array<double, 3> v)
{
    return gravity<Axis>(v[0], v[1], v[2]); 
}

/// Gas pressure - polytropic equation of state.
double pressure(std::vector<double> const& state)
{
    return KAPPA * std::pow(rho(state), GAMMA);
}

double speed_of_sound(std::vector<double> const& state)
{
    OCTOPUS_ASSERT(rho(state) > 0.0);
    OCTOPUS_ASSERT(pressure(state) >= 0.0);
    return std::sqrt(GAMMA * pressure(state) / rho(state)); 
}

///////////////////////////////////////////////////////////////////////////////
// Kernels.
struct initialize : octopus::trivial_serialization
{
    void operator()(octopus::octree_server& U) const
    {
        using std::pow;
        using std::sqrt;

        double const ei0 = 1.0;
        double const tau0 = pow(ei0, 1.0 / GAMMA);
        double const rho1 = 1.0e-10;
        double const ei1 = 1.0e-10;
        double const tau1 = pow(ei1, 1.0 / GAMMA);
    
        double const eps = 0.4;
        double const R_outer = 1.0747e-4;
        double const R_inner = R_outer*(1.0 - eps)/(1.0 + eps);

        double const h = sqrt(2.0*G*M_C*R_inner*R_outer/(R_inner + R_outer));

        double const C = 0.5*pow(h/R_inner, 2) - G*M_C/R_inner;    
  
        boost::uint64_t const gnx = octopus::config().grid_node_length;
 
        for (boost::uint64_t i = 0; i < gnx; ++i)
        {
            for (boost::uint64_t j = 0; j < gnx; ++j)
            {
                for (boost::uint64_t k = 0; k < gnx; ++k)
                {
                    double const x = U.xc(i);
                    double const y = U.yc(j);
                    double const z = U.zc(k);
  
                    // Cylindrical R.  
                    double const r = sqrt(pow(x, 2) + pow(y, 2));
    
                    if ((R_inner <= r) && (R_outer >= r))
                    {
                        double const z_max =
                            sqrt(pow(G*M_C/(0.5*pow(h/r, 2) - C), 2) - r*r);

                        if (z <= z_max)
                        {
                            double const rho_here =
                                  (0.5/KAPPA)
                                * (C + G*M_C/sqrt(r*r + z*z) - 0.5*pow(h/r, 2));

                            rho(U(i, j, k))          = rho_here;
                            momentum_x(U(i, j, k))   = -y*rho_here*h/pow(r, 2);
                            momentum_y(U(i, j, k))   = x*rho_here*h/pow(r, 2);
                            total_energy(U(i, j, k)) = ei0;
                            tau(U(i, j, k))          = tau0;
                        }

                        else
                        {
                            rho(U(i, j, k))          = rho1;
                            momentum_x(U(i, j, k))   = 0.0; 
                            momentum_y(U(i, j, k))   = 0.0;
                            total_energy(U(i, j, k)) = ei1;  
                            tau(U(i, j, k))          = tau1;
                        }
                    }
                    
                    else
                    {
                        rho(U(i, j, k))          = rho1;
                        momentum_x(U(i, j, k))   = 0.0; 
                        momentum_x(U(i, j, k))   = 0.0;
                        total_energy(U(i, j, k)) = ei1;  
                        tau(U(i, j, k))          = tau1;
                    }

                    momentum_z(U(i, j, k)) = 0.0;
                }
            }
        }
    }
};

struct enforce_outflow : octopus::trivial_serialization
{
    void operator()(octopus::face f, boost::array<double, 3> x) const
    {
        // IMPLEMENT
    } 
};

// FIXME: This should live in <octopus/science/> and be a default.
struct reflect_z : octopus::trivial_serialization
{
    void operator()(std::vector<double>& state) const
    {
        momentum_z(state) = -momentum_z(state);
    }
};

struct max_eigenvalue : octopus::trivial_serialization
{
    double operator()(
        octopus::axis a
      , std::vector<double> const& state
      , boost::array<double, 3> const& v
        ) const
    {
        using std::abs;

        switch (a)
        {
            case octopus::x_axis:
                return abs(momentum_x(state) / rho(state))
                     + speed_of_sound(state);

            case octopus::y_axis:
                return abs(momentum_y(state) / rho(state))
                     + speed_of_sound(state);

            case octopus::z_axis:
                return abs(momentum_z(state) / rho(state))
                     + speed_of_sound(state);

            default: { OCTOPUS_ASSERT(false); break; }
        }

        return 0.0;
    }
};

struct conserved_to_primitive : octopus::trivial_serialization
{
    void operator()(
        std::vector<double>& state
      , boost::array<double, 3> const& v
        ) const
    {
        total_energy(state) -= kinetic_energy(state);
        momentum_x(state)   /= rho(state);
        momentum_y(state)   /= rho(state);
        momentum_z(state)   /= rho(state);
    }
};

struct primitive_to_conserved : octopus::trivial_serialization
{
    void operator()(
        std::vector<double>& state
      , boost::array<double, 3> const& X
        ) const
    {
        momentum_x(state)   *= rho(state);
        momentum_y(state)   *= rho(state);
        momentum_z(state)   *= rho(state);
        total_energy(state) += kinetic_energy(state);
    }
};

struct source : octopus::trivial_serialization
{
    std::vector<double> operator()(
        std::vector<double> const& state
      , boost::array<double, 3> const& v
        ) const
    {
        std::vector<double> s(octopus::science().state_size);

        momentum_x(s) = rho(state)*gravity<octopus::x_axis>(v);
        momentum_y(s) = rho(state)*gravity<octopus::y_axis>(v);
        momentum_z(s) = rho(state)*gravity<octopus::z_axis>(v);

        return s;
    }
};

struct floor_state : octopus::trivial_serialization
{
    void operator()(
        std::vector<double>& state
      , boost::array<double, 3> const& v
        ) const
    {
        rho(state) = (std::max)(rho(state), rho_floor); 

        double const internal_energy
            = total_energy(state) - kinetic_energy(state);

        if (internal_energy > 0.1 * total_energy(state))
        {
            using std::pow;
            tau(state) = pow((std::max)(internal_energy, internal_energy_floor)
                                      , 1.0 / GAMMA); 
        }
    }
};

struct flux : octopus::trivial_serialization
{
    std::vector<double> operator()(
        octopus::axis a 
      , std::vector<double>& state
      , boost::array<double, 3> const& v
        ) const
    {
        double p = pressure(state);
 
        std::vector<double> fl(state);

        switch (a)
        {
            case octopus::x_axis:
            {
                // Velocity.
                double const v = momentum_x(state) / rho(state);

                for (boost::uint64_t i = 0; i < state.size(); ++i)
                    fl[i] = state[i] * v;

                momentum_x(fl)   += p;
                total_energy(fl) += v * p;

                break;
            }

            case octopus::y_axis:
            {
                // Velocity.
                double const v = momentum_y(state) / rho(state);

                for (boost::uint64_t i = 0; i < state.size(); ++i)
                    fl[i] = state[i] * v;

                momentum_y(fl)   += p;
                total_energy(fl) += v * p;

                break;
            }

            case octopus::z_axis:
            {
                // Velocity.
                double const v = momentum_z(state) / rho(state);

                for (boost::uint64_t i = 0; i < state.size(); ++i)
                    fl[i] = state[i] * v;

                momentum_z(fl)   += p;
                total_energy(fl) += v * p;

                break;
            }

            default: { OCTOPUS_ASSERT(false); break; }
        }

        return fl;
    }
};

void octopus_define_problem(octopus::science_table& sci)
{
    sci.state_size = 6;

    sci.initialize = initialize();
    sci.enforce_outflow = enforce_outflow();
    sci.reflect_z = reflect_z();
    sci.max_eigenvalue = max_eigenvalue();
    sci.conserved_to_primitive = conserved_to_primitive(); 
    sci.primitive_to_conserved = primitive_to_conserved();
    sci.source = source();
    sci.floor = floor_state();
    sci.flux = flux();  

    sci.output = octopus::single_variable_silo_writer(0, "rho");
}

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    octopus::octree_init_data root_data;
    root_data.dx = octopus::science().initial_spacestep();
    root.create_root(hpx::find_here(), root_data);

    root.apply(octopus::science().initialize);

    root.output();

    ///////////////////////////////////////////////////////////////////////////
    // Crude, temporary stepper.

    // FIXME: Proper support for adding commandline options.     
    double dt = 0.0; 
    double temporal_domain = 0.0;

    octopus::config_reader reader;

    reader
        ("3d_torus.dt", dt, 1.0e-10)
        ("3d_torus.temporal_domain", temporal_domain, 2.0)
    ;

    double time = 0.0;

    std::cout << (boost::format("dt = %.6e\n") % dt)
              << (boost::format("stepping to %.6e\n") % temporal_domain)
              << "\n";  

    while (time <= temporal_domain)
    {
        std::cout << (boost::format("time = %.6e\n") % time);

        root.step(1.0e-10);
        root.output();

        time += dt;
    } 
    
    return 0;
}

