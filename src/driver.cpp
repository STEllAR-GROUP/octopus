////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_fwd.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime.hpp>
#include <hpx/runtime/components/stubs/runtime_support.hpp>
#include <hpx/lcos/future.hpp>
#include <hpx/lcos/future_wait.hpp>

#include <octopus/driver.hpp>
#include <octopus/engine/engine_server.hpp>

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;

using hpx::components::stubs::runtime_support;

namespace octopus { extern OCTOPUS_EXPORT int default_main(variables_map& vm); }

OCTOPUS_EXPORT int main(int argc, char** argv);

int hpx_main(variables_map& vm)
{
    int result = 0;

    {
        std::cout << "Launching Octopus AMR driver...\n"
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
    
        // Find all localities supporting Octopus.
        std::vector<hpx::id_type> localities =
            hpx::find_all_localities(type);
    
        OCTOPUS_ASSERT_MSG(!localities.empty(),
                           "no localities supporting Octopus available"); 
    
        std::cout << "Found " << localities.size() << " usable localities\n";

        // FIXME: Sadly, distributing factory doesn't support constructor args
        // yet, so we have to do this by hand.
        std::vector<hpx::future<hpx::id_type, hpx::naming::gid_type> > engines;    
        engines.reserve(localities.size());

        ///////////////////////////////////////////////////////////////////////
        // Define the problem. 
        typedef void (*define_function)(octopus::science_table&);
        typedef int (*main_function)(variables_map&);

        typedef boost::function<void(define_function)> define_deleter;
        typedef boost::function<void(main_function)> main_deleter;

        // Figure out where we are.
        boost::plugin::dll this_exe(hpx::util::get_executable_filename());

        std::pair<define_function, define_deleter> define_p = 
            this_exe.get<define_function, define_deleter>
                ("octopus_define_problem");

        std::pair<main_function, main_deleter> main_p = 
            this_exe.get<main_function, main_deleter>
                ("octopus_main");

        OCTOPUS_ASSERT_MSG(define_p.first || main_p.first,
            "either octopus_define_problem or octopus_main must be defined");

        // Initialize the science table.
        octopus::science_table sci = octopus::default_science_table(); 

        if (define_p.first)
            (*define_p.first)(sci);

        for (std::size_t i = 0; i < localities.size(); ++i)
        {
            engines.push_back(
                runtime_support::create_component_async<octopus::engine_server>
                    (localities[i], cfg, sci, localities));
        }
    
        hpx::wait(engines);   

        ///////////////////////////////////////////////////////////////////////
        // Invoke user entry point or default main.
        if (main_p.first)
            result = (*main_p.first)(vm);
        else
            result = octopus::default_main(vm);
    }

    hpx::finalize();
    return result;
}

int main(int argc, char** argv)
{
    options_description cmdline("Octopus AMR Driver");

    ///////////////////////////////////////////////////////////////////////////
    // Initialize HPX.
    return hpx::init(cmdline, argc, argv); 
}

