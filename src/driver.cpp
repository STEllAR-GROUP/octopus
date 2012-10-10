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

int hpx_main(variables_map& vm)
{
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

        // TODO: Temporary filler.
        octopus::science_table sci = octopus::science_table(); 

        for (std::size_t i = 0; i < localities.size(); ++i)
        {
            engines.push_back(
                runtime_support::create_component_async<octopus::engine_server>
                    (localities[i], cfg, sci, localities));
        }
    
        hpx::wait(engines);   
    }

    int result = octopus_main(vm);
    hpx::finalize();
    return result;
}

OCTOPUS_EXPORT int main(int argc, char** argv);

int main(int argc, char** argv)
{
    options_description cmdline("Octopus AMR Driver");

    ///////////////////////////////////////////////////////////////////////////
    // Initialize HPX.
    return hpx::init(cmdline, argc, argv); 
}

