////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

// http://www.vistrails.org/index.php/User:Tohline/Apps/PapaloizouPringleTori

// Once upon a time, I thought this would just have one include :P.

#include <octopus/driver.hpp>
#include <octopus/science.hpp>
#include <octopus/engine/engine_interface.hpp>
#include <octopus/engine/ini.hpp>
#include <octopus/io/silo.hpp>
#include <octopus/octree/octree_reduce.hpp>
#include <octopus/operators/boost_array_arithmetic.hpp>
#include <octopus/filesystem.hpp>

#include <boost/process.hpp> 

#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>

// FIXME: Move shared code from the drivers into a shared object/headers.
// FIXME: Names, Proper configuration.

// Gravitation constant.
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
double KAPPA = 1.0;

// Directory to store data in (default is set below to the current path). 
std::string data_directory = "";

// The directory containing the scripts that define the simulation
// visualization.
std::string visualization_directory =
    octopus::join_paths(OCTOPUS_CURRENT_SOURCE_DIRECTORY
                      , "visualization/tablet_interactive");

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

struct cfl_timestep : octopus::trivial_serialization
{
    double operator()(octopus::octree_server& U) const
    {
        // REVIEW: I need to initialize this to some value higher than any
        // possible dt, I think...
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

// Serializable minimum
struct minimum : octopus::trivial_serialization
{
    template <typename T>
    T const& operator()(T const& a, T const& b) const
    {
        return (std::min)(a, b);
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

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    octopus::config_reader reader;

   reader
        ("visit.data_directory", data_directory,
            boost::filesystem::current_path().string())
        ("visit.visualization_directory", visualization_directory)
    ;

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

    sci.output = octopus::single_variable_silo_writer(0, "rho"
      , octopus::join_paths(data_directory, "3d_torus.silo").c_str()
        );
}

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    octopus::octree_init_data root_data;
    root_data.dx = octopus::science().initial_spacestep();
    root.create_root(hpx::find_here(), root_data);

    root.apply(octopus::science().initialize);

    root.refine();

    root.output();

    ///////////////////////////////////////////////////////////////////////////
    // Crude, temporary stepper.

    // FIXME: Proper support for adding commandline options and INI parameters. 
    double dt = 0.0; 
    double max_dt_growth = 0.0; 
    double temporal_domain = 0.0;

    octopus::config_reader reader;

    reader
        // FIXME: Move these somewhere more generic.
        ("dt", dt, 1.0e-10)
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_domain", temporal_domain, 1.0e-6)

        ("3d_torus.kappa", KAPPA, 1.0)
    ;

    std::cout
        << (boost::format("dt               = %.6e\n") % dt)
        << (boost::format("max_dt_growth    = %.6e\n") % max_dt_growth)
        << (boost::format("kappa            = %.6e\n") % KAPPA)
        << "\n"
        << (boost::format("Stepping to %.6e...\n") % temporal_domain)
        << "\n";

    double time = 0.0;

    double dt_last = 0.01*root.reduce<double>(cfl_timestep(), minimum());

    boost::uint64_t step = 0;
    double next_output_time = output_frequency;

    // Removing old sim2 file.
    //    system("rm /home/zbyerly/SC12/octopus.sim2");

    // Now launching VisIt environment.
//    std::cout << "Calling VisItSetupEnvironment()\n";
    VisItSetupEnvironment();
    //write out .sim file that VisIt uses to connect
//    std::cout << "Calling VisItInitializeSock...()\n";
    VisItInitializeSocketAndDumpSimFile(
        "3d_torus",
        "",
        "",
//        data_directory.c_str(),
//        "/srv/scratch/wash/octopus/gcc-4.6.2-debug",
        NULL,
        NULL,
        octopus::join_paths(data_directory, "3d_torus.sim2").c_str()); 

    // Hacky stuff to launch VisIt.
    //system("/opt/visit/2.5.2/bin/visit -o /home/zbyerly/SC12/octopus.sim2");
    //popen("/opt/visit/2.5.2/bin/visit -o /home/zbyerly/SC12/octopus.sim2\n");
    // FIXME: Less hardcoded.
    std::string exec = "/opt/visit/2.5.2/bin/visit";

    std::vector<std::string> args;
    //args.push_back("-fullscreen");
    args.push_back("-cli");
    args.push_back("-o");
    args.push_back(octopus::join_paths(data_directory, "3d_torus.sim2"));

    boost::process::context ctx;
    ctx.environment = boost::process::self::get_environment(); 
    // //    ctx.stdout_behavior = boost::process::silence_stream();
//     ctx.stdout_behavior = boost::process::capture_stream();
//     ctx.stderr_behavior = boost::process::capture_stream();
    
    boost::process::child c = boost::process::launch(exec, args, ctx);
    

//     boost::process::pistream &is = c.get_stderr(); 
//     std::string line; 
//     while (std::getline(is, line)) 
//       std::cout << line << std::endl; 

//    sleep(10.0);

/*    int visitstate =*/ VisItDetectInput(0, -1);
     
//    std::cout << visitstate << "\n";

    // FIXME (wash): fprintf is a great evil. 
    if(VisItAttemptToCompleteConnection()) {
      fprintf(stderr, "VisIt connected\n");
    } else {
      fprintf(stderr, "VisIt did not connect\n");
    }

    std::string initialize_script =
        octopus::join_paths(visualization_directory, "initialize.py");
    std::string reload_script =
        octopus::join_paths(visualization_directory, "reload.py");

    OCTOPUS_ASSERT(!data_directory.empty());
    VisItExecuteCommand(std::string(
        "DATA_DIRECTORY = \"" + data_directory + "\"\n" 
      + "VISUALIZATION_DIRECTORY = \"" + visualization_directory + "\"\n"
      + "Source(\"" + initialize_script + "\")\n"
    ).c_str());

    // Visit appears to delete these.
    // FIXME (wash): We need to use custom handlers here that will ensure that
    // we clean up semi-gracefully if we get killed. This is super important. 
    hpx::set_error_handlers();

    while (time < temporal_domain)
    {
        dt = (std::min)(dt_last*max_dt_growth
                      , root.reduce<double>(cfl_timestep(), minimum()));

        OCTOPUS_ASSERT(0.0 < dt);

        std::cout << ( boost::format("STEP %06u : %.6e += %.6e\n")
                     % step % time % dt);

        root.step(dt);

        time += dt;
        ++step;
        dt_last = dt;

        // Commenting out so it will output every timestep for demo.
        // if ((time + dt) >= next_output_time)
        // {   
        //     std::cout << "OUTPUT\n";
        //     root.output();
        //     next_output_time += output_frequency; 
        // }

        std::cout << "RELOADING\n";
        root.output();
    
        VisItExecuteCommand(std::string(
            "Source(\"" + reload_script + "\")\n"
        ).c_str());
    } 

    std::cout << "DONE!\n";

    VisItExecuteCommand("Close()");
    VisItExecuteCommand("quit()");

    //std::cout << "Deleted plots..\n";

    //if(!VisItProcessEngineCommand())
    //VisItDisconnect();

    //std::cout << "disconnected";

    return 0;
}

