////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
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
#include <octopus/octree/octree_reduce.hpp>
#include <octopus/octree/octree_apply_leaf.hpp>
#include <octopus/operators/boost_array_arithmetic.hpp>
#include <octopus/math.hpp>

// FIXME: Move shared code from the drivers into a shared object/headers.
// FIXME: Names.
// FIXME: Proper configuration.

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

// Polytropic constant.
double KAPPA = 1.0;
    
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
inline double gravity(boost::array<double, 3> const& v)
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

        double const kappa0 = 1.0;

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
                    double const x = U.x_center(i);
                    double const y = U.y_center(j);
                    double const z = U.z_center(k);
  
                    // Cylindrical R.  
                    double const r = sqrt(pow(x, 2) + pow(y, 2));

                    // DEBUGGING
                    //std::cout << "r       = " << r       << "\n"
                    //          << "R_inner = " << R_inner << "\n"
                    //          << "R_outer = " << R_outer << "\n";
    
                    if ((R_inner <= r) && (R_outer >= r))
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
        ) const
    {
        boost::array<double, 3> coords;
        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        return (*this)(a, state, coords);
    }

    double operator()(
        octopus::axis a
      , std::vector<double> const& state
      , boost::array<double, 3> const& 
        ) const
    {
        using std::abs;

        switch (a)
        {
            case octopus::x_axis:
                return abs(momentum_x(state) / rho(state)) // velocity
                     + speed_of_sound(state);

            case octopus::y_axis:
                return abs(momentum_y(state) / rho(state)) // velocity
                     + speed_of_sound(state);

            case octopus::z_axis:
                return abs(momentum_z(state) / rho(state)) // velocity
                     + speed_of_sound(state);

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
                    std::vector<double> const& u = U(i, j, k);
                    double const dx = U.get_dx(); 

                    // FIXME: 0.4 shouldn't be hard coded.  
                    double const dt_here_x
                        = 0.4*dx/(max_eigenvalue()(octopus::x_axis, u));
                    double const dt_here_y
                        = 0.4*dx/(max_eigenvalue()(octopus::y_axis, u));
                    double const dt_here_z
                        = 0.4*dx/(max_eigenvalue()(octopus::z_axis, u));
  
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
                                        , octopus::minimum_functor());
    }
};

// IMPLEMENT: Post prediction.
struct cfl_predict_timestep
{
  private:
    double max_dt_growth_;
    double fudge_factor_;

  public:
    cfl_predict_timestep() {}

    cfl_predict_timestep(
        double max_dt_growth
      , double fudge_factor
        )
      : max_dt_growth_(max_dt_growth)
      , fudge_factor_(fudge_factor)
    {}

    /// Returns the tuple (timestep N + 1 size, timestep N+gap size)
    octopus::timestep_prediction operator()(
        octopus::octree_server& root
        ) const
    {
        OCTOPUS_ASSERT(0 < max_dt_growth_);
        OCTOPUS_ASSERT(0 < fudge_factor_);

        boost::uint64_t const gap = octopus::config().temporal_prediction_gap;

        OCTOPUS_ASSERT(0 == root.get_level());

        double next_dt = (std::min)(
            // get_dt may block
            root.get_dt() * max_dt_growth_
          , root.reduce<double>(cfl_treewise_compute_timestep()
                              , octopus::minimum_functor())
        );

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
        std::vector<double>& u
      , boost::array<double, 3> const& v
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
        std::vector<double>& u
      , boost::array<double, 3> const& X
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
    std::vector<double> operator()(
        std::vector<double> const& u
      , boost::array<double, 3> const& coords
        ) const
    {
        std::vector<double> s(octopus::science().state_size);

        momentum_x(s) = rho(u)*gravity<octopus::x_axis>(coords);
        momentum_y(s) = rho(u)*gravity<octopus::y_axis>(coords);
        momentum_z(s) = rho(u)*gravity<octopus::z_axis>(coords);

        return s;
    }
};

struct floor_state : octopus::trivial_serialization
{
    void operator()(
        std::vector<double>& u
      , boost::array<double, 3> const& coords 
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

struct flux : octopus::trivial_serialization
{
    std::vector<double> operator()(
        octopus::axis a 
      , std::vector<double>& u
      , boost::array<double, 3> const& coords
        ) const
    {
        double p = pressure(u);
 
        std::vector<double> fl(u);

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
    bool refine(std::vector<double> const& state)
    {
        if (rho(state) > min_refine_rho)
            return true;
        else
            return false;
    }

    /// If this returns true for all regions in a point, that region will be
    /// unrefined.
    bool unrefine(std::vector<double> const& state)
    {
        // Unused currently.
        return false;
    }
};

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    double max_dt_growth = 0.0; 
    double temporal_prediction_limiter = 0.0; 

    octopus::config_reader reader("octopus.3d_torus");

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("kappa", KAPPA, 1.0)
    ;

    std::cout
        << "[octopus.3d_torus]\n"
        << ( boost::format("max_dt_growth               = %.6e\n")
           % max_dt_growth)
        << ( boost::format("temporal_prediction_limiter = %.6e\n")
           % temporal_prediction_limiter)
        << ( boost::format("kappa                       = %.6e\n")
           % KAPPA)
        << "\n";

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

    sci.initial_timestep = cfl_initial_timestep();
    sci.predict_timestep = cfl_predict_timestep
        (max_dt_growth, temporal_prediction_limiter);

    sci.refine_policy = refine_by_density();

    sci.output = octopus::single_variable_silo_writer(0, "rho");
}

struct stepper : octopus::trivial_serialization
{
    void operator()(octopus::octree_server& root) const
    {
        root.apply(octopus::science().initialize);
    
        root.refine();
    
        root.output_initial();
    
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
    
        // FIXME: Proper support for adding commandline options.     
        double max_dt_growth = 0.0; 
        double temporal_prediction_limiter = 0.0; 
    
        octopus::config_reader reader("octopus.3d_torus");
    
        reader
            ("max_dt_growth", max_dt_growth, 1.25)
            ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
            ("kappa", KAPPA, 1.0)
        ;
    
        root.post_dt(root.apply_leaf(octopus::science().initial_timestep));
        double next_output_time = octopus::config().output_frequency;
    
        while (root.get_time() < octopus::config().temporal_domain)
        {
            std::cout << ( boost::format("STEP %06u : TIME %.6e += %.6e\n")
                         % root.get_step() % root.get_time() % root.get_dt());
    
            root.step();
    
            if (root.get_time() >= next_output_time)
            {   
                std::cout << "OUTPUT\n";
                root.output();
                next_output_time += octopus::config().output_frequency; 
            }
    
            // IMPLEMENT: Futurize w/ continutation.
            octopus::timestep_prediction prediction
                = root.apply_leaf(octopus::science().predict_timestep);
    
            OCTOPUS_ASSERT(0.0 < prediction.next_dt);
            OCTOPUS_ASSERT(0.0 < prediction.future_dt);
    
            root.set_time(root.get_time() + root.get_dt());
            root.set_step(root.get_step() + 1);

            root.post_dt(prediction.next_dt);

            // Update kappa.
            // FIXME: Distributed.
            reader
                ("kappa", KAPPA, 1.0)
            ;
        } 
    }
};

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    octopus::octree_init_data root_data;
    root_data.dx = octopus::science().initial_spacestep();
    root.create_root(hpx::find_here(), root_data);

    root.apply_leaf<void>(stepper());
    
    return 0;
}

