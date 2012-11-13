////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "3d_torus.hpp"

#include <hpx/lcos/wait_all.hpp>

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

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("kappa", kappa0, 1.0)
        ("rotation_direction", rotation_direction_str, "counterclockwise")
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
        << ( boost::format("kappa                       = %.6g\n")
           % kappa0)
        << ( boost::format("rotional_direction          = %s\n")
           % rotation_direction_str)
        << "\n";

    sci.state_size = 6;

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
        root.apply(octopus::science().initialize);
  
/* 
        for ( boost::uint64_t i = 1
            ; i <= octopus::config().max_refinement_level 
            ; ++i)
        {
            root.refine(i);
            root.apply(octopus::science().initialize);
            root.child_to_parent_injection(0);
            std::cout << "refined " << i << "\n";
        }
*/

        root.refine(/*octopus::config().max_refinement_level*/);
        root.child_to_parent_injection(0);

        root.output_initial();
    
        //std::cout << "Initial state prepared\n";
    
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
    
        root.post_dt(root.apply_leaf(octopus::science().initial_timestep));
        double next_output_time = octopus::config().output_frequency;
    
        while (root.get_time() < octopus::config().temporal_domain)
        {
            char const* fmt = "STEP %06u : TIME %.6e += %.6e : KAPPA %.6g\n";

            std::cout <<
                ( boost::format(fmt)
                % root.get_step() % root.get_time() % root.get_dt()
                % kappa(root.get_step()) 
                );
    
            root.step();
    
            // Update kappa.
            hpx::wait_all(octopus::call_everywhere
                (update_kappa(root.get_step() - 1)));

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

