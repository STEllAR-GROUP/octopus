////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_9BA6055C_E7A9_4A16_8A24_B8B410AA1A14)
#define OCTOPUS_9BA6055C_E7A9_4A16_8A24_B8B410AA1A14

// http://www.vistrails.org/index.php/User:Tohline/Apps/PapaloizouPringleTori

#include <octopus/state.hpp>
#include <octopus/driver.hpp>
#include <octopus/science.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/engine/ini.hpp>
#include <octopus/octree/octree_reduce.hpp>
#include <octopus/octree/octree_apply_leaf.hpp>
#include <octopus/math.hpp>
#include <octopus/global_variable.hpp>
#include <octopus/io/multi_writer.hpp>
#include <octopus/io/fstream.hpp>

#if defined(OCTOPUS_HAVE_SILO)
    #include <octopus/io/silo.hpp>
#endif

#include <boost/format.hpp>
#include <boost/atomic.hpp>
#include <boost/math/constants/constants.hpp>

#include <hpx/include/plain_actions.hpp>

// FIXME: Move shared code from the drivers into a shared object/headers.
// FIXME: Names.
// FIXME: Proper configuration.
// FIXME: Globals are bad mkay.

enum rotational_direction
{
    rotate_clockwise,
    rotate_counterclockwise
};

enum momentum_conservation
{
    angular_momentum_conservation,
    cartesian_momentum_conservation
};

///////////////////////////////////////////////////////////////////////////////
double const initial_cfl_factor = 1.0e-2;

double const cfl_factor = 0.4;

/// Gravitation constant.
double const G = 1.0; 

/// Mass of the central object.
double const M_C = 1.0;

/// Outer radius of the torus.
double const R_outer = 1.0;

/// Polytropic index.
double const gamma_ = 1.333; 

double const polytropic_n = 3.0;

/// Amplitude of the perturbation.
double const kick_amplitude = 1e-2;

/// Angular speed of the frame, e.g. how fast the grid rotates.
OCTOPUS_GLOBAL_VARIABLE((double), omega);

/// Direction of the torus' rotation.
OCTOPUS_GLOBAL_VARIABLE((rotational_direction), rot_dir);

/// Advection scheme. 
OCTOPUS_GLOBAL_VARIABLE((momentum_conservation), mom_cons);

/// If true, a rotating frame of reference is used. 
OCTOPUS_GLOBAL_VARIABLE((bool), rotating_grid);

/// Polytropic constant.
OCTOPUS_GLOBAL_VARIABLE((double), kappa);

/// The ratio of the inner radius to the outer radius, e.g. the thinnness of the
/// torus. 
OCTOPUS_GLOBAL_VARIABLE((double), X_in);

/// Mode of the perturbation.
OCTOPUS_GLOBAL_VARIABLE((boost::uint64_t), kick_mode);

///////////////////////////////////////////////////////////////////////////////
/// Mass density
double&       rho(octopus::state& u)       { return u[0]; }
double const& rho(octopus::state const& u) { return u[0]; }

/// Momentum density (X-axis)
double&       momentum_x(octopus::state& u)       { return u[1]; }
double const& momentum_x(octopus::state const& u) { return u[1]; }

/// Momentum density (Y-axis)
double&       momentum_y(octopus::state& u)       { return u[2]; }
double const& momentum_y(octopus::state const& u) { return u[2]; }

/// Momentum density (Z-axis)
double&       momentum_z(octopus::state& u)       { return u[3]; }
double const& momentum_z(octopus::state const& u) { return u[3]; }

/// Total energy of the gas 
double&       total_energy(octopus::state& u)       { return u[4]; }
double const& total_energy(octopus::state const& u) { return u[4]; }

/// Entropy tracer
double&       tau(octopus::state& u)       { return u[5]; }
double const& tau(octopus::state const& u) { return u[5]; }

enum { radial_momentum_idx = 6 };

double&       angular_momentum(octopus::state& u)       { return u[7]; }
double const& angular_momentum(octopus::state const& u) { return u[7]; }

double radius(octopus::array<double, 3> const& v)
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1]);
}

double radius(double x, double y)
{
    return std::sqrt(x*x + y*y);
}

double radial_momentum(
    octopus::state const& u
  , octopus::array<double, 3> const& v
    )
{
    switch (mom_cons)
    {
        case angular_momentum_conservation:
            return u[radial_momentum_idx]; 
        case cartesian_momentum_conservation:
        {
            double const R = radius(v);
            return momentum_x(u)*v[0]/R - momentum_y(u)*v[1]/R;
        }
        default: break;
    }

    OCTOPUS_ASSERT(false);
    return 0.0;
}

double tangential_momentum(
    octopus::state const& u
  , octopus::array<double, 3> const& v
    )
{
    double const R = radius(v);
    
    switch (mom_cons)
    {
        case angular_momentum_conservation:
            return angular_momentum(u)/R - rho(u)*R*omega; 
        case cartesian_momentum_conservation:
            return momentum_y(u)*v[0]/R
                 - momentum_x(u)*v[1]/R
                 - rho(u)*R*omega;
        default: break;
    }

    OCTOPUS_ASSERT(false);
    return 0.0;
}

template <octopus::axis Axis>
double gravity(octopus::array<double, 3> const& v)
{
    double const x = v[0];
    double const y = v[1];
    double const z = v[2];

    double const r = std::sqrt(x*x + y*y + z*z);
    double const F = -G*M_C/(r*r);

    switch (Axis)
    {
        case octopus::x_axis:
            return F*(x/r);
        case octopus::y_axis:
            return F*(y/r);
        case octopus::z_axis:
            return F*(z/r);
        default: break;
    }

    OCTOPUS_ASSERT(false);
    return 0.0; 
}

double radial_gravity(octopus::array<double, 3> const& v)
{
    double const x = v[0];
    double const y = v[1];
    double const z = v[2];

    double const r = std::sqrt(x*x + y*y + z*z);
    double const R = radius(v);
    double const F = -G*M_C/(r*r);

    return F*(R/r);
}

template <octopus::axis Axis>
double velocity(
    octopus::state const& u
  , octopus::array<double, 3> const& v
    )
{
    double const x = v[0];
    double const y = v[1];

    double const R = radius(v);
    double const st = tangential_momentum(u, v);
    double const sr = radial_momentum(u, v);
 
    switch (Axis)
    {
        case octopus::x_axis:
            return (x*sr-y*st)/R/rho(u);
        case octopus::y_axis:
            return (y*sr+x*st)/R/rho(u);
        case octopus::z_axis:
            return momentum_z(u)/rho(u);
        default: break;
    }

    OCTOPUS_ASSERT(false);
    return 0.0; 
}

double kinetic_energy(
    octopus::state const& u
  , octopus::array<double, 3> const& v
    )
{
    double const sr = radial_momentum(u, v);
    double const st = tangential_momentum(u, v);
    return 0.5 * (sr*sr + st*st + momentum_z(u)*momentum_z(u)) / rho(u); 
}

/// Gas pressure - polytropic equation of state.
double pressure(octopus::state const& u)
{
    return kappa * std::pow(rho(u), gamma_);
}

double speed_of_sound(octopus::state const& u)
{
    OCTOPUS_ASSERT(rho(u) > 0.0);
    OCTOPUS_ASSERT(pressure(u) >= 0.0);
    return std::sqrt(gamma_ * pressure(u) / rho(u)); 
}

// Omega at the radius with the highest pressure.
double omega_R_0()
{
    double const j_H = std::sqrt(2.0*X_in/(1.0+X_in));
    double const j_here = j_H*std::sqrt(G*M_C*R_outer);
    return (G*M_C)*(G*M_C)/(j_here*j_here*j_here);
}

void initialize_omega()
{
    if (rotating_grid)
        omega = omega_R_0();
    else
        omega = 0.0; 
}

double orbital_period()
{
    double const pi = boost::math::constants::pi<double>();
    return 2.0*pi/omega_R_0(); 
}

double z_max(double R)
{
    using std::pow;
    using std::sqrt;

    double const X = R/R_outer;
    double const j_H = sqrt(2.0*X_in/(1.0+X_in));
    double const C = 1.0/(1.0+X_in);
    double const tmp = pow(C+0.5*(j_H/X)*(j_H/X), -2.0) - X*X;

    if (tmp <= 0.0)
        return 0.0; 

    return R_outer*sqrt(tmp);
}

double rho_max() 
{
    using std::pow;
    using std::sqrt;

    double const n = polytropic_n;

    double const R_0 = R_outer*2.0*X_in/(1.0+X_in);        
    double const j_H = sqrt(2.0*X_in/(1.0+X_in));
    double const X_max = R_0/R_outer;
    double const C = 1.0/(1.0+X_in);
    double const H_max = 1.0/sqrt(X_max*X_max) - 0.5*(j_H/X_max)*(j_H/X_max)-C;

    return pow(H_max/((n+1)*kappa), n);
}

double density_floor()
{
    return 1e-10 * rho_max();
}

///////////////////////////////////////////////////////////////////////////////
// Kernels.
struct initialize : octopus::trivial_serialization 
{ // {{{
    void operator()(octopus::octree_server& U) const
    {
        using std::pow;
        using std::sqrt;
        using std::atan2;
        using std::cos;

        double const ei0 = 1.0;
        double const tau0 = pow(ei0, 1.0 / gamma_);
        double const rho1 = density_floor();
        double const ei1 = density_floor();
        double const tau1 = pow(ei1, 1.0 / gamma_);
    
        double const R_inner = X_in * R_outer;

        double const C = 1.0/(1.0+X_in);

        double const j_H = sqrt(2.0*X_in/(1.0+X_in));
        // Conversion to "real" units (REVIEW: What does this mean?)
        double const j_here = j_H*sqrt(G*M_C*R_outer); 
  
        boost::uint64_t const gnx = octopus::config().grid_node_length;
 
        for (boost::uint64_t i = 0; i < gnx; ++i)
        {
            for (boost::uint64_t j = 0; j < gnx; ++j)
            {
                for (boost::uint64_t k = 0; k < gnx; ++k)
                {
                    double const x_here = U.x_center(i);
                    double const y_here = U.y_center(j);
                    // REVIEW: Why do we do std::abs() here?
                    double const z_here = std::fabs(U.z_center(k));
  
                    double const R = radius(x_here, y_here);

                    // DEBUGGING
                    //std::cout << "r       = " << R       << "\n"
                    //          << "R_inner = " << R_inner << "\n"
                    //          << "R_outer = " << R_outer << "\n";
    
                    if ((R_inner <= R) && (R_outer >= R))
                    {
                        // DEBUGGING
                        //std::cout << "z     = " << z_here << "\n"
                        //          << "z_max = " << z_max  << "\n";

                        if (z_here <= z_max(R))
                        {
                            double const X = R/R_outer;

                            double const z = z_here/R_outer;
                            double const H_here = 1.0/sqrt(X*X+z*z)
                                                - 0.5*(j_H/X)*(j_H/X) - C;

                            double const n = polytropic_n;
                            double rho_here = (pow(H_here/((n+1)*kappa), n))
                                            * (G*M_C/R_outer);

                            OCTOPUS_ASSERT(rho_here > 0.0);

                            double const theta_here = atan2(y_here, x_here);

                            if (kick_mode != 0)
                                rho_here *= 1.0 + ( kick_amplitude
                                                  * cos(kick_mode*theta_here));

                            OCTOPUS_ASSERT(rho_here > 0.0);

                            rho(U(i, j, k)) = rho_here;

                            double& mom_x = momentum_x(U(i, j, k));
                            double& mom_y = momentum_y(U(i, j, k));

                            switch (rot_dir)
                            {
                                case rotate_counterclockwise:
                                {
                                    mom_x = (-y_here)*rho_here*j_here/pow(R, 2);
                                    mom_y = (+x_here)*rho_here*j_here/pow(R, 2);
                                    break;
                                }

                                case rotate_clockwise:
                                {
                                    mom_x = (+y_here)*rho_here*j_here/pow(R, 2);
                                    mom_y = (-x_here)*rho_here*j_here/pow(R, 2);
                                    break;
                                }

                                default: OCTOPUS_ASSERT(false);
                            }

                            total_energy(U(i, j, k))     = ei0;
                            tau(U(i, j, k))              = tau0;
                            angular_momentum(U(i, j, k)) = j_here*rho_here;
                        }

                        else
                        {
                            rho(U(i, j, k))              = rho1;
                            momentum_x(U(i, j, k))       = 0.0; 
                            momentum_y(U(i, j, k))       = 0.0;
                            total_energy(U(i, j, k))     = ei1;  
                            tau(U(i, j, k))              = tau1;
                            angular_momentum(U(i, j, k)) = 0.0;
                        }
                    }
                    
                    else
                    {
                        rho(U(i, j, k))              = rho1;
                        momentum_x(U(i, j, k))       = 0.0; 
                        momentum_y(U(i, j, k))       = 0.0;
                        total_energy(U(i, j, k))     = ei1;  
                        tau(U(i, j, k))              = tau1;
                        angular_momentum(U(i, j, k)) = 0.0;
                    }

                    // DEBUGGING
                    //std::cout << "(" << x_here
                    //          << ", " << y_here
                    //          << ", " << z_here << ") == "
                    //          << rho(U(i, j, k)) << "\n";

                    momentum_z(U(i, j, k))          = 0.0; 
                    U(i, j, k)[radial_momentum_idx] = 0.0;

//                    rho(U(i, j, k)) = (std::max)(rho(U(i, j, k))
//                                               , density_floor()); 
                }
            }
        }
    }
}; // }}}

struct enforce_outflow : octopus::trivial_serialization
{
    void operator()(
        octopus::octree_server& U
      , octopus::state& u
      , octopus::array<double, 3> const& loc
      , octopus::face f
        ) const
    {
        std::cout << "ENFORCE OUTFLOW: (" << loc[0]
                                  << ", " << loc[1]
                                  << ", " << loc[2]
                                  << ") " << f;
        switch (f)
        {
            case octopus::XU:
            {
                if (velocity<octopus::x_axis>(u, loc) > 0.0)
                {
                    std::cout << " OUTFLOW";
//                    total_energy(u) -= 0.5*momentum_x(u)*momentum_x(u)/rho(u);
                    momentum_x(u) = 0.0;

                    //double const vy = velocity<octopus::y_axis>(u, loc);
                    double const R = radius(loc);
//                    angular_momentum(u) = loc[0]*velocity<octopus::y_axis>(u, loc)*rho(u);
//                    u[radial_momentum_idx] = loc[1]*velocity<octopus::y_axis>(u, loc)*rho(u)/R;
                }
                break;
            }
            case octopus::XL:
            {
                if (velocity<octopus::x_axis>(u, loc) < 0.0)
                {
                    std::cout << " OUTFLOW";
//                    total_energy(u) -= 0.5*momentum_x(u)*momentum_x(u)/rho(u);
                    momentum_x(u) = 0.0;

                    //double const vy = velocity<octopus::y_axis>(u, loc);
                    double const R = radius(loc);
//                    angular_momentum(u) = loc[0]*velocity<octopus::y_axis>(u, loc)*rho(u);
//                    u[radial_momentum_idx] = loc[1]*velocity<octopus::y_axis>(u, loc)*rho(u)/R;
                }
                break;
            }

            case octopus::YU:
            {
                if (velocity<octopus::y_axis>(u, loc) > 0.0)
                {
                    std::cout << " OUTFLOW";
//                    total_energy(u) -= 0.5*momentum_y(u)*momentum_y(u)/rho(u);
                    momentum_y(u) = 0.0;

                    //double const vx = velocity<octopus::x_axis>(u, loc);
                    double const R = radius(loc);
//                    angular_momentum(u) = -loc[1]*velocity<octopus::x_axis>(u, loc)*rho(u);
//                    u[radial_momentum_idx] = loc[0]*velocity<octopus::x_axis>(u, loc)*rho(u)/R;
                }
                break;
            }
            case octopus::YL:
            {
                if (velocity<octopus::y_axis>(u, loc) < 0.0)
                {
                    std::cout << " OUTFLOW";
//                    total_energy(u) -= 0.5*momentum_y(u)*momentum_y(u)/rho(u);
                    momentum_y(u) = 0.0;

                    //double const vx = velocity<octopus::x_axis>(u, loc);
                    double const R = radius(loc);
//                    angular_momentum(u) = -loc[1]*velocity<octopus::x_axis>(u, loc)*rho(u);
//                    u[radial_momentum_idx] = loc[0]*velocity<octopus::x_axis>(u, loc)*rho(u)/R;
                }
                break;
            }

            case octopus::ZU:
            {
                if (momentum_z(u) > 0.0)
                {
                    std::cout << " OUTFLOW";
//                    total_energy(u) -= 0.5*momentum_z(u)*momentum_z(u)/rho(u);
                    momentum_z(u) = 0.0;
                }
                break;
            }
            case octopus::ZL:
            {
                if (momentum_z(u) < 0.0)
                {
                    std::cout << " OUTFLOW";
//                    total_energy(u) -= 0.5*momentum_z(u)*momentum_z(u)/rho(u);
                    momentum_z(u) = 0.0;
                }
                break;
            }

            default: OCTOPUS_ASSERT(false); break;
        }

        std::cout << "\n";
    } 
};

struct enforce_lower_limits : octopus::trivial_serialization
{
    void operator()(
        octopus::state& u
      , octopus::array<double, 3> const& v 
        ) const
    {
        rho(u) = (std::max)(rho(u), density_floor()); 

        double const internal_energy = total_energy(u) - kinetic_energy(u, v);

        if (internal_energy > 0.1 * total_energy(u))
        {
            tau(u) = std::pow((std::max)(internal_energy, density_floor())
                                       , 1.0 / gamma_); 
        }

        // Floor everything in the center of the grid.
        double const R_inner = X_in * R_outer;

        if (radius(v) < 0.5 * R_inner)
        {
            // REVIEW: Why don't we floor tau and energy? 
            rho(u)                 = density_floor();
            momentum_x(u)          = 0.0;
            momentum_y(u)          = 0.0;
            momentum_z(u)          = 0.0;
            angular_momentum(u)    = 0.0;
            u[radial_momentum_idx] = 0.0;
        }

        // If we're not conserving angular momentum, define angular momentum 
        // in terms of x and y momentum. 
        // FIXME: Not sure this belongs in this particular function, or perhaps
        // this hook should be renamed to be something more generic than
        // enforce_lower_limits.
        else if (mom_cons != angular_momentum_conservation)
            angular_momentum(u) = momentum_y(u)*v[0] - momentum_x(u)*v[1];
    }
};

struct reflect_z : octopus::trivial_serialization
{
    void operator()(octopus::state& s) const
    {
        momentum_z(s) = -momentum_z(s);
    }
};

struct max_eigenvalue : octopus::trivial_serialization
{
    double operator()(
        octopus::octree_server& U
      , octopus::state const& u
      , octopus::array<double, 3> const& v 
      , octopus::axis a
        ) const
    {
        switch (a)
        {
            case octopus::x_axis:
                return std::fabs(velocity<octopus::x_axis>(u, v)) 
                     + speed_of_sound(u);

            case octopus::y_axis:
                return std::fabs(velocity<octopus::y_axis>(u, v)) 
                     + speed_of_sound(u);

            case octopus::z_axis:
                return std::fabs(velocity<octopus::z_axis>(u, v)) 
                     + speed_of_sound(u);

            default: { OCTOPUS_ASSERT(false); break; }
        }

        return 0.0;
    }
};

// FIXME: Refactor reconstruction harness (nearly identical code is used in
// octree_server for computing fluxes).
struct cfl_treewise_compute_dt : octopus::trivial_serialization
{
    cfl_treewise_compute_dt() {} 

    double compute_x_dt(octopus::octree_server& U) const
    { // {{{ 
        double dt_inv = 0.0; 

        boost::uint64_t const bw = octopus::science().ghost_zone_width;
        boost::uint64_t const gnx = octopus::config().grid_node_length;
    
        std::vector<octopus::state> q0(gnx, octopus::state());
        std::vector<octopus::state> ql(gnx, octopus::state());
        std::vector<octopus::state> qr(gnx, octopus::state());
    
        for (boost::uint64_t k = bw; k < (gnx - bw); ++k)
            for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
            {
                for (boost::uint64_t i = 0; i < gnx; ++i)
                {
                    q0[i] = U(i, j, k);
        
                    octopus::array<double, 3> loc = U.center_coords(i, j, k);
        
                    octopus::science().conserved_to_primitive(q0[i], loc);
                }
        
                octopus::science().reconstruct(q0, ql, qr);
        
                for (boost::uint64_t i = bw; i < gnx - bw + 1; ++i)
                {
                    octopus::array<double, 3> loc = U.x_face_coords(i, j, k);
        
                    octopus::science().primitive_to_conserved(ql[i], loc);
                    octopus::science().primitive_to_conserved(qr[i], loc);
       
                    double const l_dt_inv =
                        (max_eigenvalue()(U, ql[i], loc, octopus::x_axis));
                    double const r_dt_inv = 
                        (max_eigenvalue()(U, qr[i], loc, octopus::x_axis));

                    dt_inv = octopus::maximum(dt_inv, l_dt_inv, r_dt_inv);
                    OCTOPUS_ASSERT(0.0 < dt_inv);
                }
            }

        return cfl_factor * (1.0 / (dt_inv / (U.get_dx())));
    } // }}}
    
    double compute_y_dt(octopus::octree_server& U) const
    { // {{{ 
        double dt_inv = 0.0; 

        boost::uint64_t const bw = octopus::science().ghost_zone_width;
        boost::uint64_t const gnx = octopus::config().grid_node_length;
    
        std::vector<octopus::state> q0(gnx, octopus::state());
        std::vector<octopus::state> ql(gnx, octopus::state());
        std::vector<octopus::state> qr(gnx, octopus::state());
    
        for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
            for (boost::uint64_t k = bw; k < (gnx - bw); ++k)
            {
                for (boost::uint64_t j = 0; j < gnx; ++j)
                {
                    q0[j] = U(i, j, k);
        
                    octopus::array<double, 3> loc = U.center_coords(i, j, k);
        
                    octopus::science().conserved_to_primitive(q0[j], loc);
                }
        
                octopus::science().reconstruct(q0, ql, qr);
        
                for (boost::uint64_t j = bw; j < gnx - bw + 1; ++j)
                {
                    octopus::array<double, 3> loc = U.y_face_coords(i, j, k);
        
                    octopus::science().primitive_to_conserved(ql[j], loc);
                    octopus::science().primitive_to_conserved(qr[j], loc);
        
                    double const l_dt_inv =
                        (max_eigenvalue()(U, ql[j], loc, octopus::y_axis));
                    double const r_dt_inv = 
                        (max_eigenvalue()(U, qr[j], loc, octopus::y_axis));

                    dt_inv = octopus::maximum(dt_inv, l_dt_inv, r_dt_inv);
                    OCTOPUS_ASSERT(0.0 < dt_inv);
                }
            }

        return cfl_factor * (1.0 / (dt_inv / (U.get_dx())));
    } // }}}
    
    double compute_z_dt(octopus::octree_server& U) const
    { // {{{ 
        double dt_inv = 0.0; 

        boost::uint64_t const bw = octopus::science().ghost_zone_width;
        boost::uint64_t const gnx = octopus::config().grid_node_length;
    
        std::vector<octopus::state> q0(gnx, octopus::state());
        std::vector<octopus::state> ql(gnx, octopus::state());
        std::vector<octopus::state> qr(gnx, octopus::state());
    
        for (boost::uint64_t i = bw; i < (gnx - bw); ++i)
            for (boost::uint64_t j = bw; j < (gnx - bw); ++j)
            {
                for (boost::uint64_t k = 0; k < gnx; ++k)
                {
                    q0[k] = U(i, j, k);
        
                    octopus::array<double, 3> loc = U.center_coords(i, j, k);
    
                    octopus::science().conserved_to_primitive(q0[k], loc);
                }
        
                octopus::science().reconstruct(q0, ql, qr);
        
                for (boost::uint64_t k = bw; k < gnx - bw + 1; ++k)
                {
                    octopus::array<double, 3> loc = U.z_face_coords(i, j, k);
        
                    octopus::science().primitive_to_conserved(ql[k], loc);
                    octopus::science().primitive_to_conserved(qr[k], loc);
        
                    double const l_dt_inv =
                        (max_eigenvalue()(U, ql[k], loc, octopus::z_axis));
                    double const r_dt_inv =
                        (max_eigenvalue()(U, qr[k], loc, octopus::z_axis));

                    dt_inv = octopus::maximum(dt_inv, l_dt_inv, r_dt_inv);
                    OCTOPUS_ASSERT(0.0 < dt_inv);
                }
            }

        return cfl_factor * (1.0 / (dt_inv / (U.get_dx())));
    } // }}}

    // Compute maximum dt locally in parallel. 
    double operator()(octopus::octree_server& U) const
    {
        // Do two directions in other threads.
        boost::array<hpx::future<double>, 2> xy =
        { {
            hpx::async(boost::bind
                (&cfl_treewise_compute_dt::compute_x_dt, this, boost::ref(U)))
          , hpx::async(boost::bind
                (&cfl_treewise_compute_dt::compute_y_dt, this, boost::ref(U)))
        } };

        // And do one direction here.
        double dt_limit = compute_z_dt(U);
        OCTOPUS_ASSERT(0.0 < dt_limit);

        // Wait for the x and y computations.
        dt_limit = (std::min)(dt_limit, xy[0].move());
        OCTOPUS_ASSERT(0.0 < dt_limit);

        dt_limit = (std::min)(dt_limit, xy[1].move());
        OCTOPUS_ASSERT(0.0 < dt_limit);

        return dt_limit;
    }
};

struct cfl_initial_dt : octopus::trivial_serialization
{
    double operator()(octopus::octree_server& root) const
    {
        return initial_cfl_factor 
             * root.reduce<double>(cfl_treewise_compute_dt()
                                 , octopus::minimum_functor()
                                 , std::numeric_limits<double>::max());
    }
};

// IMPLEMENT: Post prediction.
struct cfl_predict_dt
{
  private:
    double max_dt_growth_;
    double fudge_factor_;

  public:
    cfl_predict_dt() : max_dt_growth_(0.0), fudge_factor_(0.0) {}

    cfl_predict_dt(
        double max_dt_growth
      , double fudge_factor
        )
      : max_dt_growth_(max_dt_growth)
      , fudge_factor_(fudge_factor)
    {}

    /// Returns the tuple (timestep N + 1 size, timestep N + gap size)
    octopus::dt_prediction operator()(
        octopus::octree_server& root
        ) const
    {
        OCTOPUS_ASSERT(0 < max_dt_growth_);
        OCTOPUS_ASSERT(0 < fudge_factor_);

        OCTOPUS_ASSERT(0 == root.get_level());

        double next_dt = root.reduce<double>(cfl_treewise_compute_dt()
                                           , octopus::minimum_functor()
                                           , std::numeric_limits<double>::max());

        return octopus::dt_prediction(next_dt, fudge_factor_ * next_dt); 
    }

    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & max_dt_growth_;
        ar & fudge_factor_;
    }
};

// Primitive variables are mass, velocity and pressure. Conserved variables
// are mass, momentum and energy. 

struct conserved_to_primitive : octopus::trivial_serialization
{
    void operator()(
        octopus::state& u
      , octopus::array<double, 3> const& v
        ) const
    {
        double const R = radius(v);

        total_energy(u)        -= kinetic_energy(u, v);
        momentum_x(u)          /= rho(u);
        momentum_y(u)          /= rho(u);
        momentum_z(u)          /= rho(u);
        u[radial_momentum_idx] /= rho(u);
        angular_momentum(u)    /= rho(u) * R;
        angular_momentum(u)    -= omega * R;
    }
};

struct primitive_to_conserved : octopus::trivial_serialization
{
    void operator()(
        octopus::state& u
      , octopus::array<double, 3> const& v
        ) const
    {
        double const R = radius(v);

        momentum_x(u)          *= rho(u);
        momentum_y(u)          *= rho(u);
        momentum_z(u)          *= rho(u);
        u[radial_momentum_idx] *= rho(u);
        angular_momentum(u)    += omega * R;
        angular_momentum(u)    *= rho(u) * R;
        total_energy(u)        += kinetic_energy(u, v);
    }
};

struct source : octopus::trivial_serialization
{
    octopus::state operator()(
        octopus::octree_server& U 
      , octopus::state const& u
      , octopus::array<double, 3> const& v
        ) const
    {
        octopus::state s;

/*
        if (  octopus::compare_real(v[0], -1.26562, 1e-5)
           && octopus::compare_real(v[1], -1.45312, 1e-5)
           && octopus::compare_real(v[2], 0.046875, 1e-5))
        {
            std::stringstream ss;
            ss << ( boost::format("SOURCE U (%g, %g, %g):")
                  % v[0] % v[1] % v[2]); 
            for (boost::uint64_t i = 0; i < u.size(); ++i)
                ss << ( boost::format(" %.16x")
                      % octopus::hex_real(u[i]));
            ss << "\n";
//            ss << ( boost::format("Z GRAVITY: %.17e\n")
//                  % gravity<octopus::z_axis>(v));
            std::cout << ss.str();
        }

        return s;
*/

        double const R = radius(v);
        double const p = pressure(u);
        double const lz = angular_momentum(u);

        // This is the radial momentum source term (independent of gravity).
        s[radial_momentum_idx] += (p + std::pow(lz/R, 2)/rho(u))/R;

        // Add half of the Coriolis force for rotating cartesian momentum.
        momentum_x(s) += momentum_y(u)*omega;
        momentum_y(s) -= momentum_x(u)*omega; 

        // There won't be any gravity within a certain radius of the center.
        if (R < 0.5*(X_in*R_outer))
            return s;        

        momentum_x(s) += rho(u)*gravity<octopus::x_axis>(v);
        momentum_y(s) += rho(u)*gravity<octopus::y_axis>(v);
        momentum_z(s) += rho(u)*gravity<octopus::z_axis>(v);

        s[radial_momentum_idx] += rho(u)*radial_gravity(v);

        return s;
    }
};

// REVIEW: For rot_dir to work properly, do we need to make some changes here?
struct flux : octopus::trivial_serialization
{
    octopus::state operator()(
        octopus::octree_server& U
      , octopus::state& u
      , octopus::array<double, 3> const& loc
      , octopus::array<boost::uint64_t, 3> const& idx
      , octopus::axis a 
        ) const
    {
        double const p = pressure(u);
        double const R = radius(loc);
 
        octopus::state fl;

        switch (a)
        {
            case octopus::x_axis:
            {
                // Velocity.
                double const v = velocity<octopus::x_axis>(u, loc);

                for (boost::uint64_t i = 0; i < u.size(); ++i)
                    fl[i] = u[i] * v;

                momentum_x(fl)   += p;
                total_energy(fl) += v * p;

                fl[radial_momentum_idx] += loc[0] * p / R;
                angular_momentum(fl)    -= loc[1] * p;

/*
                if (  octopus::compare_real(U.x_center(idx[0]), -1.26562, 1e-5)
                   && octopus::compare_real(U.y_center(idx[1]), -1.45312, 1e-5)
                   && octopus::compare_real(U.z_center(idx[2]), 0.046875, 1e-5))
                {
                    std::stringstream ss;
                    ss << ( boost::format("X FLUX U (%g, %g, %g):")
                          % U.x_center(idx[0])
                          % U.y_center(idx[1])
                          % U.z_center(idx[2]));
                    for (boost::uint64_t i = 0; i < fl.size(); ++i)
                        ss << ( boost::format(" %.16x")
                              % octopus::hex_real(u[i]));
                    ss << "\n";
                    ss << ( boost::format("X FLUX (%g, %g, %g):")
                          % U.x_center(idx[0])
                          % U.y_center(idx[1])
                          % U.z_center(idx[2]));
                    for (boost::uint64_t i = 0; i < fl.size(); ++i)
                        ss << ( boost::format(" %.16x")
                              % octopus::hex_real(fl[i]));
                    ss << "\n";
                    ss << (boost::format("RADIUS: %.17e\n") % R);
                    ss << (boost::format("PRESSURE: %.17e\n") % p);
                    ss << (boost::format("VELOCITY: %.17e\n") % v);
                    std::cout << ss.str();
                }
*/

                break;
            }

            case octopus::y_axis:
            {
                double const v = velocity<octopus::y_axis>(u, loc);

                for (boost::uint64_t i = 0; i < u.size(); ++i)
                    fl[i] = u[i] * v;

                momentum_y(fl)   += p;
                total_energy(fl) += v * p;

                fl[radial_momentum_idx] += loc[1] * p / R;
                angular_momentum(fl)    += loc[0] * p;

                break;
            }

            case octopus::z_axis:
            {
                double const v = velocity<octopus::z_axis>(u, loc);

                for (boost::uint64_t i = 0; i < u.size(); ++i)
                    fl[i] = u[i] * v;

                momentum_z(fl)   += p;
                total_energy(fl) += v * p;

                break;
            }

            default: { OCTOPUS_ASSERT(false); break; }
        }

        return fl;
    }
};

struct refine_by_geometry
  : octopus::elementwise_refinement_criteria_base<refine_by_geometry>
{
    /// Returns true if we should refine the region that contains this point.
    bool refine(
        octopus::octree_server& U
      , octopus::state const& u
      , octopus::array<double, 3> loc
        )
    {
        using std::sqrt;
        using std::abs;

        double const dx = U.get_dx();

        double const xU = loc[0]+0.5*dx; // x upper
        double const xL = loc[0]-0.5*dx; // x lower
        double const yU = loc[1]+0.5*dx; // y upper
        double const yL = loc[1]-0.5*dx; // y lower

        double const radius0 = radius(loc[0], loc[1]); 
        double const radius1 = radius(xU, yU);
        double const radius2 = radius(xU, yL);
        double const radius3 = radius(xL, yU);
        double const radius4 = radius(xL, yL);

        bool condition0 = (
            std::fabs(loc[2]+0.5*dx) < z_max(radius0) ||
            std::fabs(loc[2]-0.5*dx) < z_max(radius0)
            );

        bool condition1 = (
            ((radius1 < R_outer) && (radius1 > X_in*R_outer)) ||  
            ((radius2 < R_outer) && (radius2 > X_in*R_outer)) ||  
            ((radius3 < R_outer) && (radius3 > X_in*R_outer)) ||  
            ((radius4 < R_outer) && (radius4 > X_in*R_outer))
            );

        return condition0 && condition1;
    }

    /// If this returns true for all regions in a point, that region will be
    /// unrefined.
    bool unrefine(
        octopus::octree_server& U
      , octopus::state const& u
      , octopus::array<double, 3> loc
        )
    {
        // Unused currently.
        return false;
    }

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        typedef elementwise_refinement_criteria_base<refine_by_geometry>
            base_type;
        ar & hpx::util::base_object_nonvirt<base_type>(*this);
    }
};

struct output_equatorial_plane 
{
    struct slicer 
    {
      private:
        std::ofstream* ofs_;

      public:
        slicer() : ofs_(0) {}

        slicer(std::ofstream& ofs) : ofs_(&ofs) {}

        void operator()(
            octopus::octree_server& U
          , octopus::state& u
          , octopus::array<double, 3>& loc
            ) const
        {
            (*ofs_) << ( boost::format("%g %g %g %i")
                       % loc[0]
                       % loc[1]
                       % loc[2]
                       % U.get_level());

            for (boost::uint64_t i = 0; i < u.size(); ++i)
            {
                (*ofs_) << ( boost::format(" %016x")
                           % octopus::hex_real(u[i]));
            }

            (*ofs_) << "\n";
        }

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int)
        {
            OCTOPUS_ASSERT(false);
        };
    };

  private:
    octopus::axis axis_;
  public:
    output_equatorial_plane() : axis_() {}

    output_equatorial_plane(octopus::axis a) : axis_(a) {}

    void operator()(
        octopus::octree_server& U
      , std::ofstream& ofs
        ) const
    {
        U.slice_leaf(slicer(ofs), axis_, 1e-7);
    }

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & axis_;
    };
};

#endif // OCTOPUS_9BA6055C_E7A9_4A16_8A24_B8B410AA1A14

