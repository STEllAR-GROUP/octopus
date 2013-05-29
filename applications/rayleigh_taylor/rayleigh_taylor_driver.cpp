////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

// http://arxiv.org/abs/1005.0114

#include <octopus/driver.hpp>
#include <octopus/science.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/engine/ini.hpp>
#include <octopus/io/silo.hpp>
#include <octopus/octree/octree_reduce.hpp>
#include <octopus/octree/octree_apply_leaf.hpp>
#include <octopus/math.hpp>

// FIXME: Move shared code from the drivers into a shared object/headers.
// FIXME: Names.
// FIXME: Proper configuration.

// Minimum value that rho is allowed to be.
double const rho_floor = 1.0e-20;

// Minimum value that internal energy is allowed to be. 
double const internal_energy_floor = 1.0e-20;  

// Minimum refinement density.
double const min_refine_rho = 1.0e-6;

// Polytropic index.
double const GAMMA = 1.4; // EULER_GAMMA

// Polytropic constant.
double KAPPA = 1.0;

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

    // TODO: This should be tunable at runtime.
    double const g = -1.0;    

    switch (Axis)
    {
        case octopus::x_axis:
            return 0.0;
        case octopus::y_axis:
            return g;
        case octopus::z_axis:
            return 0.0;       
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
// TODO: Make equation of state switch-able at runtime.
double pressure(octopus::state const& s)
{
    OCTOPUS_ASSERT(total_energy(s) >= 0.0);
    OCTOPUS_ASSERT(kinetic_energy(s) >= 0.0);
    OCTOPUS_ASSERT(tau(s) >= 0.0);

    //Polytropic equation of state.
    //return KAPPA * std::pow(rho(state), GAMMA);

    //Ideal gas equation of state.
    double const internal_energy
        = total_energy(s) - kinetic_energy(s);
    
    if (internal_energy > 0.001 * total_energy(s))
        return (GAMMA-1.0)*internal_energy;
    else
        return (GAMMA-1.0)*(std::pow(tau(s),GAMMA));        

}

double speed_of_sound(octopus::state const& s)
{
    OCTOPUS_ASSERT(rho(s) > 0.0);
    OCTOPUS_ASSERT(pressure(s) >= 0.0);
    return std::sqrt(GAMMA * pressure(s) / rho(s)); 
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

        // "High" side.
        double const rho0 = 2.0;
        // "Low" side.
        double const rho1 = 1.0;
    
        // Initial position of the discontinuity.
        double const initial_position = 0.0; 
        
        boost::uint64_t const gnx = octopus::config().grid_node_length;

        double const y_min = octopus::config().spatial_domain;        
        double const x_min = y_min;
        
        std::cout << "y_min =" << y_min << "\n";

        double const p_base = 5.0;        

        // TODO: This should be tunable at runtime. Also, this value
        // might be wrong.
        double const g = -1.0;         
        double const pi = M_PI;
        double const big_a = 0.01;                    
        double const h = 0.005;

        for (boost::uint64_t i = 0; i < gnx; ++i)
        {
            for (boost::uint64_t j = 0; j < gnx; ++j)
            {
                for (boost::uint64_t k = 0; k < gnx; ++k)
                {
                    double const x = U.x_center(i);
                    double const y = U.y_center(j);
                    double const z = U.z_center(k);

                    bool high = true;

                    if (initial_position <= x)
                    {                        
                        double ei_here = (1.0/GAMMA)*(p_base+rho0*g*(x-x_min));
                        total_energy(U(i, j, k)) = ei_here;
                        tau(U(i, j, k))          = pow(ei_here, 1.0 / GAMMA);
                    }
                    else
                    {
                        double ei_here = (1.0/GAMMA)*(p_base+rho0*g*(x-x_min)+
                                         rho1*g*(x-initial_position));                
                        total_energy(U(i, j, k)) = ei_here;
                        tau(U(i, j, k))          = pow(ei_here, 1.0 / GAMMA);
                    }

                    double const phi_here = (big_a/2.0)*(
                        std::cos(pi*(x-x_min)/x_min)+
                        std::cos(pi*(3.0*x_min-x)/x_min)+
                        y_min);

                    rho(U(i,j,k)) = rho1 + 0.5*(rho0-rho1)*(
                        1.0+std::tanh((y-phi_here)/h));

                    momentum_x(U(i, j, k)) = 0.0; 
                    momentum_y(U(i, j, k)) = 0.0; 
                    momentum_z(U(i, j, k)) = 0.0; 

                    rho(U(i, j, k)) = (std::max)(rho(U(i, j, k)), rho_floor); 
                }
            }
        }
    }
};

struct enforce_outflow : octopus::trivial_serialization
{
    void operator()(octopus::face f, octopus::array<double, 3> x) const
    {
        // IMPLEMENT
    } 
};

// FIXME: This should live in <octopus/science/> and be a default.
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
                     + speed_of_sound(s);

            case octopus::y_axis:
                return abs(momentum_y(s) / rho(s)) // velocity
                     + speed_of_sound(s);

            case octopus::z_axis:
                return abs(momentum_z(s) / rho(s)) // velocity
                     + speed_of_sound(s);

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
        octopus::octree_server& U
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

struct flux : octopus::trivial_serialization
{
    octopus::state operator()(
        octopus::octree_server& U
      , octopus::axis a 
      , octopus::state& u
      , octopus::array<double, 3> const& coords
        ) const
    {
        double p = pressure(u);
 
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
        if (rho(s) > min_refine_rho)
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

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    double max_dt_growth = 0.0; 
    double temporal_prediction_limiter = 0.0; 

    octopus::config_reader reader("octopus.rayleigh_taylor");

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("kappa", KAPPA, 1.0)
    ;
   

    std::cout
        << "[octopus.rayleigh_taylor]\n"
        << ( boost::format("max_dt_growth               = %lf\n")
           % max_dt_growth)
        << ( boost::format("temporal_prediction_limiter = %i\n")
           % temporal_prediction_limiter)
        << ( boost::format("kappa                       = %lf\n")
           % KAPPA)
        << "\n";

    sci.initialize = initialize();
    sci.enforce_outflow = enforce_outflow();
    sci.reflect_z = reflect_z();
    sci.max_eigenvalue = max_eigenvalue();
    sci.conserved_to_primitive = conserved_to_primitive(); 
    sci.primitive_to_conserved = primitive_to_conserved();
    sci.source = source();
    sci.enforce_limits = enforce_lower_limits();
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
        for ( std::size_t i = 0
            ; i < octopus::config().levels_of_refinement
            ; ++i)
        {
            std::cout << "Refining level " << i << "\n";

            root.apply(octopus::science().initialize);
            root.refine();
            root.child_to_parent_state_injection(0);

            std::cout << "Refined level " << i << "\n";
        }

        root.output("U_L%06u_initial.silo");
    
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
    
        // FIXME: Proper support for adding commandline options.     
        double max_dt_growth = 0.0; 
        double temporal_prediction_limiter = 0.0; 
 
        octopus::config_reader reader("octopus.rayleigh_taylor");
    
        reader
            ("max_dt_growth", max_dt_growth, 1.25)
            ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
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

