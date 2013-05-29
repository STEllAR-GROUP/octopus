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

    std::string rotational_direction_str = "";

    octopus::config_reader reader("octopus.3d_torus");

    double kappa0 = 0.0;
    double eps = 0.0; 
    double R_outer = 0.0; 

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("rotational_direction", rotational_direction_str, "counterclockwise")
        ("kappa", kappa0, 1.0)
        ("eps", eps, 0.4)
        ("outer_radius", R_outer, 1.0747e-4)
    ;

    kappa_buffer.store(kappa0);
    KAPPA.resize(octopus::config().temporal_prediction_gap, kappa0);
 
    if (rotational_direction_str == "clockwise")
        rotation = rotate_clockwise;
    else if (rotational_direction_str == "counterclockwise")
        rotation = rotate_counterclockwise;
    else
        OCTOPUS_ASSERT_MSG(false, "invalid rotational direction");

    std::cout
        << "[octopus.3d_torus]\n"
        << ( boost::format("max_dt_growth               = %lf\n")
           % max_dt_growth)
        << ( boost::format("temporal_prediction_limiter = %i\n")
           % temporal_prediction_limiter)
        << ( boost::format("rotational_direction        = %s\n")
           % rotational_direction_str)
        << ( boost::format("kappa                       = %.6g\n")
           % kappa0)
        << ( boost::format("epsilon                     = %.6g\n")
           % eps)
        << ( boost::format("outer_radius                = %.6g\n")
           % R_outer)
        << "\n";

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

struct stepper : octopus::trivial_serialization
{
    void operator()(octopus::octree_server& root) const
    {
        for ( std::size_t i = 0
            ; i < octopus::config().levels_of_refinement
            ; ++i)
        {
            std::cout << "REFINING LEVEL " << i << "\n";

            root.apply(octopus::science().initialize);
            root.refine();
            root.child_to_parent_state_injection(0);

            std::cout << "REFINED LEVEL " << i << "\n";
        }

        root.output("U_L%06u_initial.silo");
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

