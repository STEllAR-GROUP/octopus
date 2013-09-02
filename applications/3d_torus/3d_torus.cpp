////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "3d_torus.hpp"

#include <fenv.h>

#include <hpx/util/high_resolution_timer.hpp>

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    double max_dt_growth = 0.0; 
    double temporal_prediction_limiter = 0.0; 

    std::string rot_dir_str = "";
    std::string mom_cons_str = "";

    octopus::config_reader reader("octopus.3d_torus");

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("rotational_direction", rot_dir_str, "counterclockwise")
        ("momentum_conservation", mom_cons_str, "angular")
        ("rotating_grid", rotating_grid, true)
        ("kappa", kappa, 1.0)
        ("X_in", X_in, 0.5)
        ("kick_mode", kick_mode, 0) 
    ;

    if (rot_dir_str == "clockwise")
        rot_dir = rotate_clockwise;
    else if (rot_dir_str == "counterclockwise")
        rot_dir = rotate_counterclockwise;
    else
        OCTOPUS_ASSERT_MSG(false, "invalid rotational direction");

    if (mom_cons_str == "angular")
        mom_cons = angular_momentum_conservation;
    else if (mom_cons_str == "cartesian")
        mom_cons = cartesian_momentum_conservation;
    else
        OCTOPUS_ASSERT_MSG(false, "invalid momentum conservation");

    std::cout
        << "[octopus.3d_torus]\n"
        << ( boost::format("max_dt_growth                 = %lf\n")
           % max_dt_growth)
        << ( boost::format("temporal_prediction_limiter   = %i\n")
           % temporal_prediction_limiter)
        << ( boost::format("rotational_direction          = %s\n")
           % rot_dir_str)
        << ( boost::format("momentum_conservation         = %s\n")
           % mom_cons_str)
        << ( boost::format("rotating_grid                 = %i\n")
           % rotating_grid.get())
        << ( boost::format("kappa                         = %.6g\n")
           % kappa.get())
        << ( boost::format("X_in                          = %.6g\n")
           % X_in.get())
        << ( boost::format("kick_mode                     = %i\n")
           % kick_mode.get())
        << "\n";

    // FIXME: Move this into core code.
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);

    initialize_omega();

    std::cout << "R_0     = " << (R_outer*2.0*X_in/(1.0+X_in)) << "\n";
    std::cout << "R_inner = " << (X_in*R_outer) << "\n";
    std::cout << "rho_max = " << rho_max() << "\n";
    std::cout << "omega   = " << omega.get() << "\n";
    std::cout << "period  = " << orbital_period() << "\n\n";

    sci.initialize = initialize();
    sci.enforce_outflow = enforce_outflow();
    sci.reflect_z = reflect_z();
    sci.max_eigenvalue = max_eigenvalue();
    sci.conserved_to_primitive = conserved_to_primitive(); 
    sci.primitive_to_conserved = primitive_to_conserved();
    sci.source = source();
    sci.enforce_limits = enforce_lower_limits();
    sci.flux = flux();  

    sci.initial_dt = cfl_initial_dt();
    sci.predict_dt = cfl_predict_dt(max_dt_growth, temporal_prediction_limiter);

    sci.refine_policy = refine_by_geometry();
    sci.distribute = slice_distribution();

/*
    octopus::multi_writer mw;

    #if defined(OCTOPUS_HAVE_SILO)
        mw.add_writer(octopus::single_variable_silo_writer(0, "rho"));
    #endif

    mw.add_writer(octopus::fstream_writer(
        output_equatorial_plane(octopus::x_axis)
      , "slice_x_L%06u_S%06u.dat"));

    mw.add_writer(octopus::fstream_writer(
        output_equatorial_plane(octopus::y_axis)
      , "slice_y_L%06u_S%06u.dat"));

    mw.add_writer(octopus::fstream_writer(
        output_equatorial_plane(octopus::z_axis)
      , "slice_z_L%06u_S%06u.dat"));

    sci.output = mw;
*/
/*
    #if defined(OCTOPUS_HAVE_SILO)
        sci.output = octopus::single_variable_silo_writer(0, "rho");
    #endif
*/
    sci.output = octopus::fstream_writer(
        output_equatorial_plane(octopus::z_axis)
      , "slice_z_L%06u_S%06u.dat");
}

struct stepper 
{
  private:
    double period_;

  public:
    stepper() : period_(0.0) {}

    stepper(double period) : period_(period) {}

    void operator()(octopus::octree_server& root) const
    {
        hpx::util::high_resolution_timer refine_clock;

        boost::uint64_t refine_passes = 0;

        if (octopus::config().levels_of_refinement != 0)
        {   
            if (octopus::config().levels_of_refinement > 1)
                refine_passes = octopus::config().levels_of_refinement;
            else
                refine_passes = 2;
        }

        root.apply(octopus::science().initialize);

        for (boost::uint64_t i = 0; i < refine_passes; ++i)
        {
            root.refine();
            root.apply(octopus::science().initialize);
            std::cout << "REFINEMENT PASS " << (i + 1)
                      << " OF " << refine_passes << ", "
                      << count_nodes(root) << " NODES" << "\n"; 
        }

        root.child_to_parent_state_injection(0);

//        #if defined(OCTOPUS_HAVE_SILO)
//            root.output(root.get_time() / period_, "U_L%06u_initial.silo");
//        #endif

//        root.communicate_ghost_zones(0);

        double refine_walltime = refine_clock.elapsed();
 
        root.output(0.0);

        std::ofstream dt_file("dt.csv");
        std::ofstream speed_file("speed.csv");
 
        //dt_file    << "step, time [orbits], dt [orbits], output & refine?\n";
        //speed_file << "step, speed [orbits/hours], output & refine?\n";
        dt_file    << "# step, time [orbits], dt [orbits], dt cfl [orbits], output?\n";
        speed_file << "# step, speed [orbits/hours], output?\n";
 
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
    
        root.post_dt(root.apply_leaf(octopus::science().initial_dt));
//        root.post_dt(initial_cfl_factor*0.001);
        double next_output_time = octopus::config().output_frequency * period_;

        hpx::reset_active_counters();

        hpx::util::high_resolution_timer global_clock;
   
        bool last_step = false;

        while (!last_step)
        {
            hpx::util::high_resolution_timer local_clock;

            boost::uint64_t const this_step = root.get_step();
            double const this_dt = root.get_dt();
            double const this_time = root.get_time();

            if (  ((this_time + this_dt) / period_)
               >= octopus::config().temporal_domain)
                last_step = true;

            root.step();
   
            bool output_and_refine = false;

            if (root.get_time() >= next_output_time)
            {   
                output_and_refine = true;

                root.output(root.get_time() / period_);
                next_output_time +=
                    (octopus::config().output_frequency * period_); 

                //root.refine();
            }
  
            // IMPLEMENT: Futurize w/ continutation.
            octopus::dt_prediction prediction
                = root.apply_leaf(octopus::science().predict_dt);
 
//            octopus::dt_prediction prediction(0.001, 0.001);
   
            OCTOPUS_ASSERT(0.0 < prediction.next_dt);
            OCTOPUS_ASSERT(0.0 < prediction.future_dt);

            double next_dt = std::min(prediction.next_dt, root.get_dt() * 1.25);

            if (!last_step && (  ((root.get_time() + next_dt) / period_)
                              >= octopus::config().temporal_domain))
            {
                double t = octopus::config().temporal_domain * period_ - root.get_time();
                OCTOPUS_ASSERT(t <= next_dt);
                next_dt = t;
//                std::cout << (boost::format("TIME: %.17e\n") % root.get_time());
            }

            root.post_dt(next_dt);

            ///////////////////////////////////////////////////////////////////
            // I/O of stats
            char const* fmt = "STEP %06u : ORBITS %.7g %|34t| += %.7g "
                              "%|52t|: SPEED %.7g %|76t| [orbits/hour] ";

            double const speed =
                ((this_dt / period_) / (local_clock.elapsed() / 3600));

            std::cout <<
                ( boost::format(fmt)
                % this_step
                % (this_time / period_)
                % (this_dt / period_)
                % speed 
                );
 
            //if (output_and_refine)
            //    std::cout << ": OUTPUT & REFINE";
            if (output_and_refine)
                std::cout << ": OUTPUT";

            std::cout << "\n";
 
            // Record timestep size.
            dt_file << ( boost::format("%i %e %e %e %i\n")
                       % this_step 
                       % (this_time / period_) 
                       % (this_dt / period_) 
                       % (prediction.next_dt / period_)
                       % output_and_refine); 

            // Record speed. 
            speed_file << ( boost::format("%e %e %i\n")
                          % this_step 
                          % speed 
                          % output_and_refine); 
        }

        double solve_walltime = global_clock.elapsed();

        std::cout << "\n"
                  << "REFINE WALLTIME " << refine_walltime << " [seconds]\n" 
                  << "SOLVE WALLTIME  " << solve_walltime << " [seconds]\n" 
                  << "TOTAL WALLTIME  "
                  << (refine_walltime + solve_walltime)
                  << " [seconds]\n"; 
    }

    template <typename Archive>
    void serialize(Archive& ar, unsigned int)
    {
        ar & period_;
    }
};

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    octopus::octree_init_data root_data;
    // FIXME: create_root or root_data should do this.
    root_data.dx = octopus::science().initial_dx();
    root.create_root(hpx::find_here(), root_data);

    root.apply_leaf<void>(stepper(orbital_period()));
    
    return 0;
}

