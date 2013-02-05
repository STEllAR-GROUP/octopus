////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "3d_torus.hpp"

#include <hpx/util/high_resolution_timer.hpp>

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    double max_dt_growth = 0.0; 
    double temporal_prediction_limiter = 0.0; 

    std::string rotation_direction_str = "";

    octopus::config_reader reader("octopus.3d_torus");

    double kappa0 = 0.0;
    double eps = 0.0; 
    double R_outer = 0.0; 

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("rotation_direction", rotation_direction_str, "counterclockwise")
        ("kappa", kappa0, 1.0)
        ("epsilon", eps, 0.4)
        ("outer_radius", R_outer, 1.0747e-4)
    ;

    kappa_buffer.store(kappa0);
    KAPPA.resize(octopus::config().temporal_prediction_gap, kappa0);
 
    if (rotation_direction_str == "clockwise")
        rotation = rotate_clockwise;
    else if (rotation_direction_str == "counterclockwise")
        rotation = rotate_counterclockwise;
    else
        OCTOPUS_ASSERT_MSG(false, "invalid rotation direction");

    std::cout
        << "[octopus.3d_torus]\n"
        << ( boost::format("max_dt_growth               = %lf\n")
           % max_dt_growth)
        << ( boost::format("temporal_prediction_limiter = %i\n")
           % temporal_prediction_limiter)
        << ( boost::format("rotional_direction          = %s\n")
           % rotation_direction_str)
        << ( boost::format("kappa                       = %.6g\n")
           % kappa0)
        << ( boost::format("epsilon                     = %.6g\n")
           % eps)
        << ( boost::format("outer_radius                = %.6g\n")
           % R_outer)
        << "\n";

    sci.state_size = 6;

    sci.initialize = initialize(eps, R_outer);
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

struct stepper 
{
  private:
    double period_;

  public:
    stepper() : period_(0.0) {}

    stepper(double period) : period_(period) {}

    void operator()(octopus::octree_server& root) const
    {
        hpx::util::high_resolution_timer global_clock;
   
        for ( std::size_t i = 0
            ; i < octopus::config().levels_of_refinement
            ; ++i)
        {
            root.apply(octopus::science().initialize);
            root.refine();

            root.child_to_parent_injection(0);
            root.prepare_compute_queues();

            std::cout << "REFINED LEVEL " << (i + 1) << "\n";
        }

        root.output(root.get_time() / period_, "U_L%06u_initial.silo");
  
        std::ofstream dt_file("dt.csv");
        std::ofstream speed_file("speed.csv");
 
        dt_file    << "step, time [orbits], dt [orbits], output & refine?\n";
        speed_file << "step, speed [orbits/hours], output & refine?\n";
 
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
    
        root.post_dt(root.apply_leaf(octopus::science().initial_timestep));
        double next_output_time = octopus::config().output_frequency * period_;

        while ((root.get_time() / period_) <= octopus::config().temporal_domain)
        {
            hpx::util::high_resolution_timer local_clock;

            boost::uint64_t const this_step = root.get_step();
            double const this_dt = root.get_dt();
            double const this_time = root.get_time();

            root.step();
   
            bool output_and_refine = false;

            if (root.get_time() >= next_output_time)
            {   
                output_and_refine = true;

                root.output(root.get_time() / period_);
                next_output_time +=
                    (octopus::config().output_frequency * period_); 

                root.refine();
            }
  
            // IMPLEMENT: Futurize w/ continutation.
            octopus::timestep_prediction prediction
                = root.apply_leaf(octopus::science().predict_timestep);
    
            OCTOPUS_ASSERT(0.0 < prediction.next_dt);
            OCTOPUS_ASSERT(0.0 < prediction.future_dt);
    
            root.post_dt(prediction.next_dt);

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
 
            if (output_and_refine)
                std::cout << ": OUTPUT & REFINE";

            std::cout << "\n";
 
            // Record timestep size.
            dt_file << this_step << ", "
                    << (this_time / period_) << ", "
                    << (this_dt / period_) << ", "
                    << output_and_refine << std::endl; 

            // Record speed. 
            speed_file << this_step << ", "
                       << speed << ", "
                       << output_and_refine << std::endl;
        }

        std::cout << "\n"
                  << "ELAPSED WALLTIME "
                  << global_clock.elapsed()
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
    root_data.dx = octopus::science().initial_spacestep();
    root.create_root(hpx::find_here(), root_data);

    double eps = 0.0; 
    double R_outer = 0.0; 

    octopus::config_reader reader("octopus.3d_torus");

    reader
        ("epsilon", eps, 0.4)
        ("outer_radius", R_outer, 1.0747e-4)
    ;

    root.apply_leaf<void>(stepper(orbital_period(eps, R_outer)));
    
    return 0;
}

