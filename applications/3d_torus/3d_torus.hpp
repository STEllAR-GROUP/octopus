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

#if defined(OCTOPUS_HAVE_SILO)
    #include <octopus/io/silo.hpp>
#endif

#include <boost/atomic.hpp>
#include <boost/math/constants/constants.hpp>

#include <hpx/include/plain_actions.hpp>

// FIXME: Move shared code from the drivers into a shared object/headers.
// FIXME: Names.
// FIXME: Proper configuration.
// FIXME: Globals are bad mkay.

// Gravitation constant.
double const G = 1.0; // BIG_G

// Mass of the central object.
double const M_C = 2e-2;

// Minimum value that rho is allowed to be.
double const rho_floor = 1.0e-20;

// Minimum value that internal energy is allowed to be. 
double const internal_energy_floor = 1.0e-20;  

// Minimum refinement density.
double const min_refine_rho = 1.0e-6;

// Polytropic index.
double const GAMMA = 2.0; // EULER_GAMMA

enum rotational_direction
{
    rotate_clockwise,
    rotate_counterclockwise
};

///////////////////////////////////////////////////////////////////////////////
// Polytropic constant.
std::vector<double> KAPPA = std::vector<double>();

boost::atomic<double> kappa_buffer(0.0);

inline double& kappa(boost::uint64_t step)
{
    return KAPPA[step % octopus::config().temporal_prediction_gap];
}

void update_kappa(double k)
{
    std::cout << "Updating kappa => " << k << "\n";
    kappa_buffer.store(k);
}
HPX_PLAIN_ACTION(update_kappa, update_kappa_action);
HPX_ACTION_HAS_CRITICAL_PRIORITY(update_kappa_action);

struct set_kappa_from_buffer
{
  private:
    boost::uint64_t step_;

  public:
    set_kappa_from_buffer() : step_(0) {}

    set_kappa_from_buffer(boost::uint64_t step) : step_(step) {}

    void operator()() const
    {
        kappa(step_) = kappa_buffer;
    }
 
    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & step_;
    }
};

struct initialize_kappa
{
  private:
    double kappa0_;

  public:
    initialize_kappa() : kappa0_(0.0) {}

    initialize_kappa(double kappa0) : kappa0_(kappa0) {}

    void operator()() const
    {
        kappa_buffer.store(kappa0_);
        KAPPA.resize(octopus::config().temporal_prediction_gap, kappa0_);
    }

    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & kappa0_;
    }
};

///////////////////////////////////////////////////////////////////////////////
/// Mass density
inline double&       rho(octopus::state& u)       { return u[0]; }
inline double const& rho(octopus::state const& u) { return u[0]; }

/// Momentum density (X-axis)
inline double&       momentum_x(octopus::state& u)       { return u[1]; }
inline double const& momentum_x(octopus::state const& u) { return u[1]; }

/// Momentum density (Y-axis)
inline double&       momentum_y(octopus::state& u)       { return u[2]; }
inline double const& momentum_y(octopus::state const& u) { return u[2]; }

/// Momentum density (Z-axis)
inline double&       momentum_z(octopus::state& u)       { return u[3]; }
inline double const& momentum_z(octopus::state const& u) { return u[3]; }

/// Total energy of the gas 
inline double&       total_energy(octopus::state& u)       { return u[4]; }
inline double const& total_energy(octopus::state const& u) { return u[4]; }

/// Entropy tracer
inline double&       tau(octopus::state& u)       { return u[5]; }
inline double const& tau(octopus::state const& u) { return u[5]; }

inline double kinetic_energy(octopus::state const& s)
{
    return 0.5 * ( momentum_x(s) * momentum_x(s)
                 + momentum_y(s) * momentum_y(s)
                 + momentum_z(s) * momentum_z(s)) / rho(s);
}

template <octopus::axis Axis>
inline double gravity(double x, double y, double z)
{
    using std::sqrt;

    double const r = sqrt(x*x + y*y + z*z);
    double const F = -G*M_C/(r*r);

    switch (Axis)
    {
        case octopus::x_axis:
            return F*(x/r);
        case octopus::y_axis:
            return F*(y/r);
        case octopus::z_axis:
        {
            double const r_cyl = sqrt(x*x + y*y); 
            return F*(z/r_cyl);
        }
        default: break;
    }

    OCTOPUS_ASSERT(false);
    return 0.0; 
}

template <octopus::axis Axis>
inline double gravity(octopus::array<double, 3> const& v)
{
    return gravity<Axis>(v[0], v[1], v[2]); 
}

/// Gas pressure - polytropic equation of state.
double pressure(octopus::state const& s, boost::uint64_t step)
{
    return kappa(step) * std::pow(rho(s), GAMMA);
}

double speed_of_sound(octopus::state const& s, boost::uint64_t step)
{
    OCTOPUS_ASSERT(rho(s) > 0.0);
    OCTOPUS_ASSERT(pressure(s, step) >= 0.0);
    return std::sqrt(GAMMA * pressure(s, step) / rho(s)); 
}

double orbital_period(double eps, double R_outer)
{
    double const R_inner = R_outer*(1.0 - eps)/(1.0 + eps);
    double const h = sqrt(2.0*G*M_C*R_inner*R_outer/(R_inner + R_outer));
    double const omega_R_max = (G*M_C)*(G*M_C)/(h*h*h);
    double const pi = boost::math::constants::pi<double>();

    return 2.0*pi/omega_R_max; 
}

///////////////////////////////////////////////////////////////////////////////
// Kernels.
struct initialize 
{ // {{{
  private:
    double eps_;
    double R_outer_;
    rotational_direction rotation_;

  public:
    initialize() : eps_(0.0), R_outer_(0.0), rotation_() {}

    initialize(
        double eps
      , double R_outer
      , rotational_direction rotation
        )
      : eps_(eps)
      , R_outer_(R_outer)
      , rotation_(rotation)
    {}

    void operator()(octopus::octree_server& U) const
    {
        using std::pow;
        using std::sqrt;

        double const kappa0 = 1.0;

        double const ei0 = 1.0;
        double const tau0 = pow(ei0, 1.0 / GAMMA);
        double const rho1 = rho_floor;
        double const ei1 = internal_energy_floor;
        double const tau1 = pow(ei1, 1.0 / GAMMA);
    
        double const R_inner = R_outer_*(1.0 - eps_)/(1.0 + eps_);

        double const h = sqrt(2.0*G*M_C*R_inner*R_outer_/(R_inner + R_outer_));

        double const C = 0.5*pow(h/R_inner, 2) - G*M_C/R_inner;    
  
        boost::uint64_t const gnx = octopus::config().grid_node_length;
 
        for (boost::uint64_t i = 0; i < gnx; ++i)
        {
            for (boost::uint64_t j = 0; j < gnx; ++j)
            {
                for (boost::uint64_t k = 0; k < gnx; ++k)
                {
                    double const x = U.x_center(i);
                    double const y = U.y_center(j);
                    double const z = U.z_center(k);
  
                    // Cylindrical R.  
                    double const r = sqrt(pow(x, 2) + pow(y, 2));

                    // DEBUGGING
                    //std::cout << "r       = " << r       << "\n"
                    //          << "R_inner = " << R_inner << "\n"
                    //          << "R_outer = " << R_outer_ << "\n";
    
                    if ((R_inner <= r) && (R_outer_ >= r))
                    {
                        double const z_max =
                            sqrt(pow(G*M_C/(0.5*pow(h/r, 2) - C), 2) - r*r);

                        // DEBUGGING
                        //std::cout << "z     = " << z     << "\n"
                        //          << "z_max = " << z_max << "\n";

                        if (z <= z_max)
                        {
                            double const rho_here =
                                  (0.5/kappa0)
                                * (C + G*M_C/sqrt(r*r + z*z) - 0.5*pow(h/r, 2));

                            rho(U(i, j, k)) = rho_here;

                            double& mom_x = momentum_x(U(i, j, k));
                            double& mom_y = momentum_y(U(i, j, k));

                            mom_x = y*rho_here*h/pow(r, 2);
                            mom_y = x*rho_here*h/pow(r, 2);

                            if (rotate_counterclockwise == rotation_)
                            {
                                mom_x = (-y)*rho_here*h/pow(r, 2);
                                mom_y = (+x)*rho_here*h/pow(r, 2);
                            }

                            else if (rotate_clockwise == rotation_)
                            {
                                mom_x = (+y)*rho_here*h/pow(r, 2);
                                mom_y = (-x)*rho_here*h/pow(r, 2);
                            }

                            else
                                OCTOPUS_ASSERT(false);

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
                        momentum_y(U(i, j, k))   = 0.0;
                        total_energy(U(i, j, k)) = ei1;  
                        tau(U(i, j, k))          = tau1;
                    }

                    // DEBUGGING
                    //std::cout << "(" << x << ", " << y << ", " << z << ") == "
                    //          << rho(U(i, j, k)) << "\n";

                    momentum_z(U(i, j, k)) = 0.0; 

                    rho(U(i, j, k)) = (std::max)(rho(U(i, j, k)), rho_floor); 
                }
            }
        }
    }

    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & eps_;
        ar & R_outer_;
    }
}; // }}}

struct enforce_outflow : octopus::trivial_serialization
{
    void operator()(octopus::face f, octopus::array<double, 3> x) const
    {
        // IMPLEMENT
    } 
};

struct enforce_lower_limits : octopus::trivial_serialization
{
    void operator()(
        octopus::state& u
      , octopus::array<double, 3> const& coords 
        ) const
    {
        rho(u) = (std::max)(rho(u), rho_floor); 

        double const internal_energy
            = total_energy(u) - kinetic_energy(u);

        if (internal_energy > 0.1 * total_energy(u))
        {
            using std::pow;
            tau(u) = pow((std::max)(internal_energy, internal_energy_floor)
                                  , 1.0 / GAMMA); 
        }
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
      , octopus::axis a
      , octopus::state const& s
        ) const
    {
        octopus::array<double, 3> coords;
        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        return (*this)(U, a, s, coords);
    }

    double operator()(
        octopus::octree_server& U
      , octopus::axis a
      , octopus::state const& s
      , octopus::array<double, 3> const& 
        ) const
    {
        using std::abs;

        switch (a)
        {
            case octopus::x_axis:
                return abs(momentum_x(s) / rho(s)) // velocity
                     + speed_of_sound(s, U.get_step());

            case octopus::y_axis:
                return abs(momentum_y(s) / rho(s)) // velocity
                     + speed_of_sound(s, U.get_step());

            case octopus::z_axis:
                return abs(momentum_z(s) / rho(s)) // velocity
                     + speed_of_sound(s, U.get_step());

            default: { OCTOPUS_ASSERT(false); break; }
        }

        return 0.0;
    }
};

struct cfl_treewise_compute_timestep : octopus::trivial_serialization
{
    double operator()(octopus::octree_server& U) const
    {
        // REVIEW (zach): I need to initialize this to some value higher than
        // any possible dt, I think...
        // REVIEW: Should this be hardcoded. 
        double dt_limit = 100.0;

        boost::uint64_t const gnx = octopus::config().grid_node_length;
        boost::uint64_t const bw = octopus::science().ghost_zone_width;

        for (boost::uint64_t i = bw; i < (gnx-bw); ++i)
        {
          for (boost::uint64_t j = bw; j < (gnx-bw); ++j)
            {
              for (boost::uint64_t k = bw; k < (gnx-bw); ++k)
                {
                    octopus::state const& u = U(i, j, k);
                    double const dx = U.get_dx(); 

                    // FIXME: 0.4 shouldn't be hard coded.  
                    double const dt_here_x
                        = 0.4*dx/(max_eigenvalue()(U, octopus::x_axis, u));
                    double const dt_here_y
                        = 0.4*dx/(max_eigenvalue()(U, octopus::y_axis, u));
                    double const dt_here_z
                        = 0.4*dx/(max_eigenvalue()(U, octopus::z_axis, u));
  
                    dt_limit = (std::min)(dt_limit, dt_here_x);
                    OCTOPUS_ASSERT(0.0 < dt_limit);
                    dt_limit = (std::min)(dt_limit, dt_here_y);
                    OCTOPUS_ASSERT(0.0 < dt_limit);
                    dt_limit = (std::min)(dt_limit, dt_here_z);
                    OCTOPUS_ASSERT(0.0 < dt_limit);
                }
            }
        }

        return dt_limit;
    }
};

struct cfl_initial_timestep : octopus::trivial_serialization
{
    double operator()(octopus::octree_server& root) const
    {
        return 0.01 * root.reduce<double>(cfl_treewise_compute_timestep()
                                        , octopus::minimum_functor()
                                        , 100.0);
    }
};

// IMPLEMENT: Post prediction.
struct cfl_predict_timestep
{
  private:
    double max_dt_growth_;
    double fudge_factor_;

  public:
    cfl_predict_timestep() : max_dt_growth_(0.0), fudge_factor_(0.0) {}

    cfl_predict_timestep(
        double max_dt_growth
      , double fudge_factor
        )
      : max_dt_growth_(max_dt_growth)
      , fudge_factor_(fudge_factor)
    {}

    /// Returns the tuple (timestep N + 1 size, timestep N + gap size)
    octopus::timestep_prediction operator()(
        octopus::octree_server& root
        ) const
    {
        OCTOPUS_ASSERT(0 < max_dt_growth_);
        OCTOPUS_ASSERT(0 < fudge_factor_);

        OCTOPUS_ASSERT(0 == root.get_level());

        double next_dt = root.reduce<double>(cfl_treewise_compute_timestep()
                                           , octopus::minimum_functor()
                                           , root.get_dt() * max_dt_growth_);

        return octopus::timestep_prediction(next_dt, fudge_factor_ * next_dt); 
    }

    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & max_dt_growth_;
        ar & fudge_factor_;
    }
};

struct conserved_to_primitive : octopus::trivial_serialization
{
    void operator()(
        octopus::state& u
      , octopus::array<double, 3> const& v
        ) const
    {
        total_energy(u) -= kinetic_energy(u);
        momentum_x(u)   /= rho(u);
        momentum_y(u)   /= rho(u);
        momentum_z(u)   /= rho(u);
    }
};

struct primitive_to_conserved : octopus::trivial_serialization
{
    void operator()(
        octopus::state& u
      , octopus::array<double, 3> const& X
        ) const
    {
        momentum_x(u)   *= rho(u);
        momentum_y(u)   *= rho(u);
        momentum_z(u)   *= rho(u);
        total_energy(u) += kinetic_energy(u);
    }
};

struct source : octopus::trivial_serialization
{
    octopus::state operator()(
        octopus::octree_server& 
      , octopus::state const& u
      , octopus::array<double, 3> const& coords
        ) const
    {
        octopus::state s;

        momentum_x(s) = rho(u)*gravity<octopus::x_axis>(coords);
        momentum_y(s) = rho(u)*gravity<octopus::y_axis>(coords);
        momentum_z(s) = rho(u)*gravity<octopus::z_axis>(coords);

        return s;
    }
};

struct flux : octopus::trivial_serialization
{
    octopus::state operator()(
        octopus::octree_server& U
      , octopus::axis a 
      , octopus::state& u
      , octopus::array<double, 3> const& coords
        ) const
    {
        double p = pressure(u, U.get_step());
 
        octopus::state fl(u);

        switch (a)
        {
            case octopus::x_axis:
            {
                // Velocity.
                double const v = momentum_x(u) / rho(u);

                for (boost::uint64_t i = 0; i < u.size(); ++i)
                    fl[i] = u[i] * v;

                momentum_x(fl)   += p;
                total_energy(fl) += v * p;

                break;
            }

            case octopus::y_axis:
            {
                // Velocity.
                double const v = momentum_y(u) / rho(u);

                for (boost::uint64_t i = 0; i < u.size(); ++i)
                    fl[i] = u[i] * v;

                momentum_y(fl)   += p;
                total_energy(fl) += v * p;

                break;
            }

            case octopus::z_axis:
            {
                // Velocity.
                double const v = momentum_z(u) / rho(u);

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

struct refine_by_density
  : octopus::elementwise_refinement_criteria_base<refine_by_density>
{
    /// Returns true if we should refine the region that contains this point.
    bool refine(octopus::state const& s)
    {
        if (rho(s) >= min_refine_rho)
            return true;
        else
            return false;
    }

    /// If this returns true for all regions in a point, that region will be
    /// unrefined.
    bool unrefine(octopus::state const& s)
    {
        // Unused currently.
        return false;
    }
};

#endif // OCTOPUS_9BA6055C_E7A9_4A16_8A24_B8B410AA1A14

