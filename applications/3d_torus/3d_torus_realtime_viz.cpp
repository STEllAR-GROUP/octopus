////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "3d_torus.hpp"

#include <octopus/filesystem.hpp>
#include <octopus/visit/visit_simulation_client.hpp>

#include <hpx/lcos/wait_all.hpp>

octopus::visit_simulation_client vsc; 

bool loop_viz;

void set_zoom(std::string const& arg)
{
    std::cout << "Setting zoom => " << arg << "\n";
    vsc.evaluate("set_zoom(" + arg + ")\n"); 
}
HPX_PLAIN_ACTION(set_zoom, set_zoom_action);
HPX_ACTION_HAS_CRITICAL_PRIORITY(set_zoom_action);

void set_angle(std::string const& arg1, std::string const& arg2)
{
    std::cout << "Setting angle => (" << arg1 << ", " << arg2 << ")\n";
    vsc.evaluate("set_angle(" + arg1 + ", " + arg2 + ")\n"); 
}
HPX_PLAIN_ACTION(set_angle, set_angle_action);
HPX_ACTION_HAS_CRITICAL_PRIORITY(set_angle_action);

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    double max_dt_growth = 0.0; 
    double temporal_prediction_limiter = 0.0; 

    std::string rotation_direction_str = "";

    octopus::config_reader torus_reader("octopus.3d_torus");

    double kappa0 = 0.0;

    torus_reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("kappa", kappa0, 1.0)
        ("rotation_direction", rotation_direction_str, "counterclockwise")
    ;
        
    octopus::config_reader visit_reader("octopus.visit");

    std::string data_directory;
    std::string visualization_directory;

    visit_reader
        ("data_directory", data_directory,
            boost::filesystem::current_path().string())
        ("visualization_directory", visualization_directory,
            octopus::join_paths(OCTOPUS_CURRENT_SOURCE_DIRECTORY
                              , "visualization"
                              , "tablet_interactive"))
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
        << ( boost::format("loop                        = %i\n")
           % loop_viz)
        << "\n";

    std::cout
        << "[octopus.visit]\n"
        << ( boost::format("data_directory          = %s\n")
           % data_directory)
        << ( boost::format("visualization_directory = %s\n")
           % visualization_directory)
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

    sci.output = octopus::single_variable_silo_writer(0, "rho"
      , octopus::join_paths(data_directory, "3d_torus.silo").c_str()
      , octopus::join_paths(data_directory, "3d_torus.silo").c_str()
        );
}

struct stepper : octopus::trivial_serialization
{
    void operator()(octopus::octree_server& root) const
    {
        octopus::config_reader visit_reader("octopus.visit");

        std::string data_directory;
        std::string visualization_directory;

        visit_reader
            ("data_directory", data_directory,
                boost::filesystem::current_path().string())
            ("visualization_directory", visualization_directory,
                octopus::join_paths(OCTOPUS_CURRENT_SOURCE_DIRECTORY
                                  , "visualization"
                                  , "tablet_interactive"))
        ;

        root.apply(octopus::science().initialize);

        root.refine();
        root.child_to_parent_injection(0);

        root.output_initial();
    
        //std::cout << "Initial state prepared\n";
    
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
    
        root.post_dt(root.apply_leaf(octopus::science().initial_timestep));
   
        vsc.create(hpx::find_here());
 
        vsc.start("3d_torus");

        std::string initialize_script =
            octopus::join_paths(visualization_directory, "initialize.py");
        std::string reload_script =
            octopus::join_paths(visualization_directory, "reload.py");

        OCTOPUS_ASSERT(!data_directory.empty());
        vsc.evaluate(
            "DATA_DIRECTORY = \"" + data_directory + "\"\n" 
          + "VISUALIZATION_DIRECTORY = \"" + visualization_directory + "\"\n"
          + "Source(\"" + initialize_script + "\")\n"
        );

        // Visit appears to delete these.
        // FIXME (wash): We need to use custom handlers here that will ensure
        // that we clean up semi-gracefully if we get killed. This is super
        // important. 
        hpx::set_error_handlers();

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
                (set_kappa_from_buffer(root.get_step() - 1)));

            root.output();

            vsc.evaluate("Source(\"" + reload_script + "\")\n");
    
            // IMPLEMENT: Futurize w/ continutation.
            octopus::timestep_prediction prediction
                = root.apply_leaf(octopus::science().predict_timestep);
    
            OCTOPUS_ASSERT(0.0 < prediction.next_dt);
            OCTOPUS_ASSERT(0.0 < prediction.future_dt);
    
            root.post_dt(prediction.next_dt);
        } 

        vsc.evaluate("Close()\n");

//        vsc.terminate();
    }
};

int octopus_main(boost::program_options::variables_map& vm)
{
    octopus::octree_client root;

    if (loop_viz)
    {
        while (true)
        {
            octopus::octree_init_data root_data;
            root_data.dx = octopus::science().initial_spacestep();
            root.create_root(hpx::find_here(), root_data);

            root.apply_leaf<void>(stepper());
        }
    }

    else
    {
        octopus::octree_init_data root_data;
        root_data.dx = octopus::science().initial_spacestep();
        root.create_root(hpx::find_here(), root_data);

        root.apply_leaf<void>(stepper());
    }
    
    return 0;
}

