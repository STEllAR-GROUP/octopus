////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//  Copyright (c) 2007-2012 Hartmut Kaiser 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_fwd.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/components/distributing_factory/distributing_factory.hpp>

#include <octopus/driver.hpp>
#include <octopus/engine/engine_server.hpp>

#include <boost/assign/std/vector.hpp>

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;

int hpx_main(variables_map& vm)
{
    {
        std::cout << "Launching Octopus AMR driver...\n\n";
    
        ///////////////////////////////////////////////////////////////////////
        // Initialize the Octopus engine.
        std::cout << "- Deploying Octopus services\n"; 
    
        hpx::components::component_type type =
            hpx::components::get_component_type<octopus::engine_server>();
    
        // Find all localities supporting Octopus.
        std::vector<hpx::id_type> localities =
            hpx::find_all_localities(type);
    
        if (localities.empty())
        {
            std::cout << "-- ERROR: No localities support Octopus\n";
            hpx::finalize();
            return 1;
        }
    
        std::cout << "-- Found " << localities.size() << " usable localities\n";
    
        typedef hpx::components::distributing_factory distributing_factory;
    
        distributing_factory factory;
        factory.create(hpx::find_here());
    
        distributing_factory::result_type result =
            factory.create_components(type, localities.size());
    
        distributing_factory::iterator_range_type parts =
            hpx::util::locality_results(result);
    
        localities.clear();
        BOOST_FOREACH(hpx::id_type const& id, parts)
        {
            localities.push_back(id);
        }
    }

    int result = octopus_main(vm);
    hpx::finalize();
    return result;
}

OCTOPUS_EXPORT int main(int argc, char** argv);

int main(int argc, char** argv)
{
    options_description cmdline("Octopus AMR Driver");

    cmdline.add_options()
        // TODO: This isn't implemented.
        //("oct:modules", "list all available Octopus modules")
        ("oct:select", value<std::vector<std::string> >()->composing(),
         "use the specified Octopus module")
        ;

    ///////////////////////////////////////////////////////////////////////////
    // Default INI
    using namespace boost::assign;
    std::vector<std::string> cfg;

    cfg += "[octopus]",
           "module_path=$[system.executable_prefix]/lib/hpx/octopus"
        ;

    ///////////////////////////////////////////////////////////////////////////
    // Initialize HPX.
    return hpx::init(cmdline, argc, argv, cfg); 
}

