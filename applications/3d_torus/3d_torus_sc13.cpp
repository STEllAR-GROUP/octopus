////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "3d_torus.hpp"

#include <fenv.h>

#include <boost/process.hpp>
#include <boost/circular_buffer.hpp>

#include <hpx/util/high_resolution_timer.hpp>
#include <hpx/runtime/threads/executors/service_executor.hpp>

#include <octopus/filesystem.hpp>

#include <boost/spirit/include/qi_rule.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_kleene.hpp>
#include <boost/spirit/include/qi_char.hpp>
#include <boost/spirit/include/qi_alternative.hpp>
#include <boost/spirit/include/qi_sequence.hpp>
#include <boost/spirit/include/qi_difference.hpp>
#include <boost/spirit/include/qi_eol.hpp>
#include <boost/spirit/include/qi_lit.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/spirit/include/support_ascii.hpp>

// Please forgive me! 
std::string gnuplot_script = ""; 
std::string buffer_directory = "";

boost::atomic<advection_scheme> new_adv_scheme(invalid_advection_scheme);

void set_advection_scheme(std::string const& arg)
{
    if (arg == "angular")
        new_adv_scheme.store(angular_advection_scheme);
    else if (arg == "cartesian")
        new_adv_scheme.store(cartesian_advection_scheme);
    else
        OCTOPUS_ASSERT_MSG(false, "invalid advection scheme");
}
HPX_PLAIN_ACTION(set_advection_scheme, set_advection_scheme_action);

void octopus_define_problem(
    boost::program_options::variables_map& vm
  , octopus::science_table& sci
    )
{
    double max_dt_growth = 0.0; 
    double temporal_prediction_limiter = 0.0; 

    std::string rot_dir_str = "";
    std::string adv_scheme_str = "";

    octopus::config_reader reader("octopus.3d_torus");

    reader
        ("max_dt_growth", max_dt_growth, 1.25)
        ("temporal_prediction_limiter", temporal_prediction_limiter, 0.5)
        ("rotational_direction", rot_dir_str, "counterclockwise")
        ("advection_scheme", adv_scheme_str, "angular")
        ("rotating_grid", rotating_grid, true)
        ("kappa", kappa, 1.0)
        ("X_in", X_in, 0.5)
        ("kick_mode", kick_mode, 0)
 
        ("sc13.gnuplot_script", gnuplot_script,
            octopus::join_paths(OCTOPUS_CURRENT_SOURCE_DIRECTORY, "sc13.gpi"))
        ("sc13.buffer_directory", buffer_directory, "/tmp/octopus_sc13_buffer")
    ;

    if (rot_dir_str == "clockwise")
        rot_dir = rotate_clockwise;
    else if (rot_dir_str == "counterclockwise")
        rot_dir = rotate_counterclockwise;
    else
        OCTOPUS_ASSERT_MSG(false, "invalid rotational direction");

    if (adv_scheme_str == "angular")
        adv_scheme = angular_advection_scheme;
    else if (adv_scheme_str == "cartesian")
        adv_scheme = cartesian_advection_scheme;
    else
        OCTOPUS_ASSERT_MSG(false, "invalid advection scheme");

    // Make sure new_adv_scheme is set to something other than invalid.
    advection_scheme expected = invalid_advection_scheme;
    new_adv_scheme.compare_exchange_strong(expected, adv_scheme);

    std::cout
        << "[octopus.3d_torus]\n"
        << ( boost::format("max_dt_growth                 = %lf\n")
           % max_dt_growth)
        << ( boost::format("temporal_prediction_limiter   = %i\n")
           % temporal_prediction_limiter)
        << ( boost::format("rotational_direction          = %s\n")
           % rot_dir_str)
        << ( boost::format("advection_scheme              = %s\n")
           % adv_scheme_str)
        << ( boost::format("rotating_grid                 = %i\n")
           % rotating_grid.get())
        << ( boost::format("kappa                         = %.6g\n")
           % kappa.get())
        << ( boost::format("X_in                          = %.6g\n")
           % X_in.get())
        << ( boost::format("kick_mode                     = %i\n")
           % kick_mode.get())
        << "\n";

    std::cout
        << "[octopus.3d_torus.sc13]\n"
        << ( boost::format("gnuplot_script                = %s\n")
           % gnuplot_script)
        << ( boost::format("buffer_directory              = %s\n")
           % buffer_directory)
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

void generate_jpeg(
    boost::uint64_t step
  , double time
  , double period
    )
{
    std::vector<std::string> args;

    const char* definitions = 
        "locality=%1%;"
        "step=%2%;"
        "time=%3%;"
        "lor=%4%;"
        "buffer_directory='%5%';"
        "suffix='L%1$06u_S%2$06u';"
        ;

    args.push_back("-e");

    // FIXME: Not sure boost.format is robust enough for this, I'd like to use
    // something that gives the parameters named variables in the gnuplot
    // script. I guess we could format them, and then stick in a couple of
    // -e statements defining them as gnuplot variables.
    args.push_back(boost::str(boost::format(definitions)
                  % hpx::get_locality_id()
                  % step
                  % (time / period)
                  % octopus::config().levels_of_refinement
                  % buffer_directory));

    args.push_back(gnuplot_script);

    boost::process::context ctx;
//    ctx.stdout_behavior = boost::process::silence_stream();
    ctx.stdout_behavior = boost::process::capture_stream();
    ctx.stderr_behavior = boost::process::capture_stream();

    std::string const exe = "/usr/bin/gnuplot";

    boost::process::child c = boost::process::launch(exe, args, ctx);

    boost::process::status s = c.wait();

    boost::process::pistream &is = c.get_stderr(); 
    std::string line; 
    while (std::getline(is, line)) 
        std::cout << line << std::endl; 

    OCTOPUS_ASSERT(!(s.exited() ? s.exit_status() : 0));
}

double relative_change(double a, double a_not)
{
    return (a - a_not) / a_not;
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
        OCTOPUS_ASSERT(gnuplot_script);
        OCTOPUS_ASSERT(buffer_directory);

        // Make sure the directory exists.
        boost::filesystem::create_directories(buffer_directory);

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

        hpx::threads::executors::io_pool_executor io_sched;

        // time, state
        typedef std::pair<double, total_state> total_state_entry;

        boost::circular_buffer<total_state_entry> cons_history(10);

        total_state_entry initial_sum(0.0, sum_state(root));

        cons_history.push_back(initial_sum);
 
        // Insert it twice so that we show up in the first few cons graphs.
        cons_history.push_back(cons_history.back());

        std::ofstream cons_file("cons.csv");

        cons_file  << "# step, time [orbits], total density, total angular momentum, \n"
                   << std::flush; 
 
        if (octopus::config().load_checkpoint)
        {
            boost::uint64_t step = 0;
            double time = 0.0;
            double dt = 0.0;

            octopus::checkpoint().read((char*) &step, sizeof(step));
            octopus::checkpoint().read((char*) &time, sizeof(time));
            octopus::checkpoint().read((char*) &dt, sizeof(dt));

            root.set_time(time, step);
            root.post_dt(dt);

            root.load();
        }

        else
        {
            root.output(0.0);
            generate_jpeg(0, 0, period_);

            std::cout << (boost::format(
                "\n"
                "TOTAL DENSITY:          %.7e\n"
                "TOTAL ANGULAR MOMENTUM: %.7e\n"
                "\n")
                % cons_history[0].second.rho
                % cons_history[0].second.angular_momentum);

                cons_file <<
                    ( boost::format("%i %e %.7e %.7e %.7e\n")
                    % 0 
                    % 0.0 
                    % cons_history[0].second.rho
                    % cons_history[0].second.angular_momentum
                    % 0.0) 
                    << std::flush;

            // Write out the "fake" cons graph for the initial state.            
            std::ofstream cons_dat_file(boost::str(boost::format(
                "cons_L%06u_S%06u.dat")
                % hpx::get_locality_id()
                % 0));
    
            cons_dat_file << (boost::format("%e %.7e\n") % 0.0 % 0.0)
                          << std::flush;
        }

        // FIXME: These aren't actually csv files anymore. 
        std::ofstream dt_file("dt.csv");
        std::ofstream speed_file("speed.csv");

        //dt_file    << "step, time [orbits], dt [orbits], output & refine?\n";
        //speed_file << "step, speed [orbits/hours], output & refine?\n";
        dt_file    << "# step, time [orbits], dt [orbits], dt cfl [orbits], output?\n"
                   << std::flush;
        speed_file << "# step, speed [orbits/hours], output?\n"
                   << std::flush;
 
        ///////////////////////////////////////////////////////////////////////
        // Crude, temporary stepper.
   
        if (octopus::config().load_checkpoint) 
        {
            // IMPLEMENT: Futurize w/ continutation.
            octopus::dt_prediction prediction
                = root.apply_leaf(octopus::science().predict_dt);
   
            OCTOPUS_ASSERT(0.0 < prediction.next_dt);
            OCTOPUS_ASSERT(0.0 < prediction.future_dt);

            double next_dt = std::min(prediction.next_dt, root.get_dt() * 1.25);

            root.post_dt(next_dt);
        }
        else
            root.post_dt(root.apply_leaf(octopus::science().initial_dt));

//        root.post_dt(initial_cfl_factor*0.001);
        double next_output_time = octopus::config().output_frequency * period_;

        hpx::reset_active_counters();

        hpx::util::high_resolution_timer global_clock;

        hpx::future<void> jpeg_future;
   
        bool last_step = false;

        while (!last_step)
        {
            // Set new advection scheme.
            advection_scheme updated = new_adv_scheme.load();

            // Do a comparison to avoid the global update if it's not necessary.
            if (updated != adv_scheme)
            {
                OCTOPUS_ASSERT(updated != invalid_advection_scheme);
                OCTOPUS_ASSERT(adv_scheme != invalid_advection_scheme);

                if (updated == angular_advection_scheme)
                    std::cout << "SWITCHING TO CARTESIAN ADVECTION SCHEME\n";  
                else 
                    std::cout << "SWITCHING TO ANGULAR ADVECTION SCHEME\n";  
                    
                adv_scheme = updated; 
            }

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
                total_state_entry ts
                    (root.get_time() / period_, sum_state(root));

                cons_file <<
                    ( boost::format("%i %e %.7e %.7e %.7e\n")
                    % this_step 
                    % (this_time / period_) 
                    % ts.second.rho
                    % ts.second.angular_momentum
                    % relative_change(
                          ts.second.angular_momentum
                        , initial_sum.second.angular_momentum))
//                        , cons_history.back().second.angular_momentum))
                    << std::flush;

                cons_history.push_back(ts);

                output_and_refine = true;

                root.output(root.get_time() / period_);

                std::ofstream cons_dat_file(boost::str(boost::format(
                    "cons_L%06u_S%06u.dat")
                    % hpx::get_locality_id()
                    % root.get_step()));

                for (boost::uint64_t i = 0; i < (cons_history.size() - 1); ++i)
                {
                    cons_dat_file <<
                        ( boost::format("%e %.7e 1e-16 -1e-16\n")
                        % cons_history[i + 1].first 
                        % relative_change(
                            cons_history[i].second.angular_momentum
                          , cons_history[i + 1].second.angular_momentum))
//                          , initial_sum.second.angular_momentum))
                        << std::flush;
                }

                if (jpeg_future.valid())
                    jpeg_future.then(boost::bind(&generate_jpeg
                      , root.get_step(), root.get_time(), period_));
                else
                    jpeg_future = hpx::async(boost::bind(&generate_jpeg
                      , root.get_step(), root.get_time(), period_));

                next_output_time +=
                    (octopus::config().output_frequency * period_); 

/*
                reset_checkpoint rc;
                hpx::wait(octopus::call_everywhere(rc));

                boost::uint64_t step = root.get_step();
                double time = root.get_time();
                double dt = root.get_dt(); 

                octopus::checkpoint().write((const char*) &step, sizeof(step));
                octopus::checkpoint().write((const char*) &time, sizeof(time));
                octopus::checkpoint().write((const char*) &dt, sizeof(dt));

                root.save();

                octopus::backup_checkpoint(".bak");
*/

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
                       % output_and_refine) 
                    << std::flush; 

            // Record speed. 
            speed_file << ( boost::format("%e %e %i\n")
                          % this_step 
                          % speed 
                          % output_and_refine)
                       << std::flush; 
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

