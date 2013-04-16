////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <hpx/hpx_fwd.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime.hpp>
#include <hpx/runtime/components/stubs/runtime_support.hpp>
#include <hpx/lcos/future.hpp>
#include <hpx/lcos/future_wait.hpp>

#include <octopus/driver.hpp>
#include <octopus/engine/engine_server.hpp>
#include <octopus/engine/engine_interface.hpp>

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;

using hpx::components::stubs::runtime_support;

// namespace octopus { extern OCTOPUS_EXPORT int default_main(variables_map& vm); }

int hpx_main(variables_map& vm)
{
    int result = 0;

    {
        // For great justice.
        std::cout << 
            "    ___          \n"
            "   (   \\  \\      \n"
            " /  \\   \\  | \\   Octopus: a scalable HPX framework for AMR\n"
            " \\__(0  0)_/ /   \n"
            "   _//||\\\\__/       Copyright (c) 2012 Bryce Adelstein-Lelbach\n"
            "  / | |\\ \\___/                         Zach Byerly\n"
            "  \\ | \\ \\__                            Dominic Marcello\n"
            "    \\_ \\_        \n"
            "\n"
            ;
         
        std::cout << "Launching AMR driver...\n"
                     "\n";

        ///////////////////////////////////////////////////////////////////////
        // Read configuration.
        octopus::config_data cfg = octopus::config_from_ini();

        std::cout << cfg << "\n"
                  "\n";
 
        ///////////////////////////////////////////////////////////////////////
        // Initialize the Octopus engine.
        hpx::components::component_type type =
            hpx::components::get_component_type<octopus::engine_server>();
    
        // Find all localities supporting Octopus (except for us).
        std::vector<hpx::id_type> localities = hpx::find_all_localities(type);

        OCTOPUS_ASSERT_MSG(!localities.empty(),
            "no localities supporting Octopus found");
    
        std::cout << "Found " << localities.size()
                  << " usable localities, deploying infrastructure...\n";

        // FIXME: Sadly, distributing factory doesn't support constructor args
        // yet, so we have to do this by hand.
        std::vector<hpx::future<hpx::id_type> > engines;    
        engines.reserve(localities.size());

        ///////////////////////////////////////////////////////////////////////
        // Define the problem. 
/*
        typedef void (*define_function)(octopus::science_table&);
        typedef int (*main_function)(variables_map&);

        typedef boost::function<void(define_function)> define_deleter;
        typedef boost::function<void(main_function)> main_deleter;

        // Figure out where we are.
        boost::plugin::dll this_exe(hpx::util::get_executable_filename());

        std::pair<define_function, define_deleter> define_p;
        std::pair<main_function, main_deleter> main_p; 

        try
        {
            main_p = this_exe.get<main_function, main_deleter>
                ("octopus_main");
        }
        catch (std::logic_error& le)
        {
            // Ignore the failure.
        }

        try
        {
            define_p = this_exe.get<define_function, define_deleter>
                ("octopus_define_problem");
        }
        catch (std::logic_error& le)
        {
            // Ignore the failure.
        }

        OCTOPUS_ASSERT_MSG(define_p.first || main_p.first,
            "either octopus_define_problem or octopus_main must be defined");
*/
        hpx::id_type const here = hpx::find_here();

        engines.emplace_back(
            runtime_support::create_component<octopus::engine_server>
                (here, cfg, octopus::default_science_table(), localities));

        ///////////////////////////////////////////////////////////////////////
        // Initialize the science table.
        std::cout << "Initializing science table...\n"
                  << "\n";

        octopus_define_problem(vm, octopus::science());

/*
        if (define_p.first)
            (*define_p.first)(sci);
*/

        ///////////////////////////////////////////////////////////////////////
        // Create an engine on every locality.
        std::cout << "Creating system components...\n";

        bool found_here = false;

        for (std::size_t i = 0; i < localities.size(); ++i)
        {
            if (localities[i] == here)
            {
                found_here = true;
                continue;
            }

            engines.emplace_back(
                runtime_support::create_component_async<octopus::engine_server>
                    (localities[i], octopus::config(),
                        octopus::science(), localities));
        }

        OCTOPUS_ASSERT(found_here);

        hpx::wait(engines);   

        ///////////////////////////////////////////////////////////////////////
        // Invoke user entry point or default main.
        std::cout << "Executing application...\n"
                     "\n";

        result = octopus_main(vm);

/*
        if (main_p.first)
            result = (*main_p.first)(vm);
        else
            result = octopus::default_main(vm);
*/
    }

    hpx::finalize();
    return result;
}

int main(int argc, char** argv)
{
    options_description cmdline("Octopus AMR Driver");

    ///////////////////////////////////////////////////////////////////////////
    // Initialize HPX.
    int r = hpx::init(cmdline, argc, argv); 

#if !defined(BOOST_MSVC)
    // We call C99 _Exit to work around problems with 3rd-party libraries using
    // atexit (HDF5 and visit).
    ::_Exit(r);
#else
    return r;
#endif
}

