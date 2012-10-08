////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//  Copyright (c) 2007-2012 Hartmut Kaiser 
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_init.hpp>
#include <hpx/runtime/naming/name.hpp>
#include <hpx/runtime/components/runtime_support.hpp>
#include <hpx/util/filesystem_compatibility.hpp>
#include <hpx/util/parse_command_line.hpp>
#include <octopus/assert.hpp>
#include <octopus/driver.hpp>
#include <octopus/engine/module.hpp>
#include <octopus/engine/engine_server.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/tokenizer.hpp>
#include <boost/plugin.hpp>

#include <set>
#include <map>
#include <iostream>

namespace fs = boost::filesystem;

namespace octopus
{

engine_server* engine_ptr = 0;

struct module
{
    std::string name_;
    boost::plugin::dll dll_;
    boost::shared_ptr<boost::plugin::plugin_factory<module_base> > factory_;
    boost::shared_ptr<module_base> handle_; 
};

/// \brief Discovers the Octopus module paths.
///
/// \returns A vector of paths. If the vector is empty, no valid paths were
///          found.
std::vector<std::string> find_module_paths()
{
    ///////////////////////////////////////////////////////////////////////////
    // Get the raw module path from the HPX INI subsystem. 
    std::string raw_path(hpx::get_config_entry("octopus.module_path",
            std::string(HPX_PREFIX) + "/lib/hpx/octopus"));

    std::cout << "-- Raw module path is '" << raw_path << "'\n";

    std::vector<std::string> mod_paths;

    ///////////////////////////////////////////////////////////////////////////
    // Split the raw path, which is delimited by : on Linux/Mac and ;
    // on Windows.
    boost::char_separator<char> sep(HPX_INI_PATH_DELIMITER);

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_type;
    tokenizer_type tok(raw_path, sep);
    tokenizer_type::iterator end = tok.end();

    for (tokenizer_type::iterator it = tok.begin(); it != end; ++it)
    {
        // If the path doesn't exist, skip it.
        if (!fs::exists(fs::path(*it)))
        {
            std::cout << "-- WARNING: Module path '" << *it
                      << "' does not exist\n";
            continue;
        }

        mod_paths.push_back(*it);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Make sure we found at least one path. 
    if (mod_paths.empty())
    {
        std::cout << "-- ERROR: Could not find any valid paths in '"
                  << raw_path << "'\n";
        return std::vector<std::string>();
    }

    ///////////////////////////////////////////////////////////////////////////
    // Report the paths that we found.
    std::cout << "-- Found module paths:\n";
    for (std::size_t i = 0; i < mod_paths.size(); ++i) 
        std::cout << "--- " << mod_paths[i] << "\n";

    return mod_paths;
}

std::set<module_role> link_modules(
    std::vector<std::string> const& requested
  , std::vector<std::string> const& mod_paths 
  , std::map<module_role, module>& mods
  , std::set<module_role>& not_fulfilled
    )
{ 
    std::set<std::string> not_found(requested.begin(), requested.end());

    for (std::size_t i = 0; i < mod_paths.size(); ++i) 
    {
        for (std::size_t j = 0; j < requested.size(); ++j)
        {
            fs::path dll_path(mod_paths[i]);
            dll_path /= fs::path(HPX_MAKE_DLL_STRING(requested[i]));

            if (!fs::exists(dll_path))
                continue;

            module mod;

            try
            {
                using boost::plugin::dll;
                using boost::plugin::plugin_factory;

                mod.name_ = requested[i];
                mod.dll_ = dll(dll_path.string(), "octopus");
                mod.factory_.reset(new plugin_factory<module_base>
                    (mod.dll_, "module"));
                mod.handle_.reset(mod.factory_->create(mod.name_));
            }
            catch (std::logic_error const& e)
            {
                std::cout << "-- WARNING: Could not dynamically link "
                          << mod.name_ << " (" << dll_path.string()
                          << "):\n--- " << e.what() << "\n"; 
                continue;
            }

            module_role r = mod.handle_->role();

            if (mods.count(r))
            {
                std::cout << "-- WARNING: Cannot use " << mod.name_
                          << " (" << dll_path.string() << ") because "
                          << mods[r].name_ << " (" << mods[r].dll_.get_name()
                          << ") already provides the " << r << " role";
                continue;
            }

            not_found.erase(requested[j]); 
            not_fulfilled.erase(r);
        }
    }

    if (!not_found.empty())
    {
        std::cout << "-- WARNING: Unable to link modules:\n";

        BOOST_FOREACH(std::string const& s, not_found)
        {
            std::cout << "--- " << s << "\n";
        }
    }

    return not_fulfilled;
}

engine_server::engine_server()
{
    OCTOPUS_ASSERT_MSG(engine_ptr == NULL,
                       "engine_ptr has already been set");
    engine_ptr = this;

    std::cout << "Starting Octopus AMR services on locality "
              << hpx::get_locality_id() << "...\n\n";
}

bool engine_server::prepare_modules()
{
    boost::program_options::variables_map vm;
    hpx::util::retrieve_commandline_arguments("Octopus AMR Driver", vm);

    std::vector<std::string> requested;
    if (vm.count("oct:select"))
        requested = vm["oct:select"].as<std::vector<std::string> >();

    ///////////////////////////////////////////////////////////////////////////
    // Find and dynamically link modules.
    std::cout << "- Preparing dynamically linked modules\n"; 
 
    std::map<module_role, module> mods;

    std::vector<std::string> mod_paths = find_module_paths();

    if (mod_paths.empty())
        // Error is logged in find_module_path().
        return false;

    std::set<module_role> not_fulfilled;
    for (std::size_t i = 0; i < last_role; ++i)
        not_fulfilled.insert(module_role(i));

    // User-specified modules.
    if (!requested.empty())
    {
        std::cout << "-- Linking user-requested modules:\n";
        for (std::size_t i = 0; i < requested.size(); ++i)
            std::cout << "--- " << requested[i] << "\n";

        link_modules(requested, mod_paths, mods, not_fulfilled);
    }

    // Defaults modules, if needed.
    if (!not_fulfilled.empty() || requested.empty())
    {
        std::vector<std::string> defaults;

        BOOST_FOREACH(module_role const& r, not_fulfilled)
        { 
            defaults.push_back(default_module(r));
        }

        std::cout << "-- Linking default modules:\n";
        for (std::size_t i = 0; i < defaults.size(); ++i)
            std::cout << "--- " << defaults[i] << "\n";

        link_modules(defaults, mod_paths, mods, not_fulfilled);
    }

    // If we still have unfilled modules, bail.   
    if (!not_fulfilled.empty())
    {
        std::cout << "-- ERROR: Some roles have not been fulfilled:\n";

        BOOST_FOREACH(module_role const& r, not_fulfilled)
        {
            std::cout << "--- " << r << "\n";
        }

        return false;
    }

    return true;
}

}

