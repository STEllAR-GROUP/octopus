////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/assert.hpp>

#include <hpx/hpx_init.hpp>
#include <hpx/include/plain_actions.hpp>
#include <hpx/async.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

///////////////////////////////////////////////////////////////////////////////
// Dummy definitions

// Merges with the running config on the server side, and enqueues a reload of 
// the configuration into the running application.
void update_live_config(hpx::util::section const& ini) { OCTOPUS_ALWAYS_ASSERT(false); }
HPX_PLAIN_ACTION(update_live_config, update_live_config_action);

// NOTE: This needs to determine whether the entry is a section or a value.
hpx::util::section get_live_config(std::string const& query)
{
    OCTOPUS_ALWAYS_ASSERT(false); 
    return hpx::util::section();
}
HPX_PLAIN_ACTION(get_live_config, get_live_config_action);

///////////////////////////////////////////////////////////////////////////////
void ini_to_json_recurse(
    std::string& json
  , hpx::util::section const& ini
    );

void ini_to_json_kernel(
    std::string& json
  , hpx::util::section const& ini
    );

std::string ini_to_json(hpx::util::section const& ini)
{
    std::string json;

    json += "{";

    ini_to_json_recurse(json, ini);

    OCTOPUS_ALWAYS_ASSERT(json.back() == ',');
    json.resize(json.size() - 1);

    json += "}";

    return json;
}

void ini_to_json_recurse(
    std::string& json
  , hpx::util::section const& ini
    )
{
    typedef hpx::util::section::section_map section_map;

    ini_to_json_kernel(json, ini);

    section_map const& sections = ini.get_sections();

    section_map::const_iterator send = sections.end();
    for (section_map::const_iterator i = sections.begin(); i != send; ++i)
        ini_to_json_recurse(json, i->second);
}

void ini_to_json_kernel(
    std::string& json
  , hpx::util::section const& ini
    )
{
    typedef hpx::util::section::entry_map entry_map;

    std::string const prefix = ini.get_full_name();

    entry_map const& entries = ini.get_entries();

    entry_map::const_iterator eend = entries.end();
    for (entry_map::const_iterator i = entries.begin(); i != eend; ++i)
    {
        std::string const expansion = ini.expand(i->second);

        json += "\"";

        if (!prefix.empty())
            json += prefix + "." + i->first;
        else
            json += i->first;

        json += "\":\"";
        json += ini.expand(i->second);
        json += "\",";
    }
}

std::vector<std::string> cgi_query_to_ini(std::string const& cgi_query)
{
    // FIXME: We should do escaping just to be safe. I don't think we need it
    // for the demo.
    std::vector<std::string> ini;
    boost::algorithm::split(ini, cgi_query,
        boost::algorithm::is_any_of("&"),
        boost::algorithm::token_compress_on);
    return ini;
} 

///////////////////////////////////////////////////////////////////////////////
void command_interpreter(hpx::id_type const& target, std::string raw_cmd)
{
    boost::algorithm::trim(raw_cmd);

    std::vector<std::string> cmd;
    boost::algorithm::split(cmd, raw_cmd,
        boost::algorithm::is_any_of(" \t\n"),
        boost::algorithm::token_compress_on);

    OCTOPUS_ALWAYS_ASSERT(!cmd.empty() && !cmd[0].empty());

    ///////////////////////////////////////////////////////////////////////////
    // Update the INI config. 
    if (cmd[0] == "update_config")
    {
        OCTOPUS_ALWAYS_ASSERT(cmd.size() == 2);

        std::vector<std::string> raw_ini = cgi_query_to_ini(cmd[1]); 

        hpx::util::section ini;
        ini.parse("runtime update", raw_ini, false);

        hpx::async<update_live_config_action>(target, ini).get();
    }

    ///////////////////////////////////////////////////////////////////////////
    // Gets the current config
    else if (cmd[0] == "get_config")
    {
        OCTOPUS_ALWAYS_ASSERT(cmd.size() == 2);

        hpx::util::section ini
            = hpx::async<get_live_config_action>(target, cmd[1]).get();

        std::cout << ini_to_json(ini) << std::flush;
    }

    else
        OCTOPUS_ALWAYS_ASSERT(false); // Unknown command.
}

///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map& vm)
{
    {
        OCTOPUS_ALWAYS_ASSERT(vm.count("target"));
        OCTOPUS_ALWAYS_ASSERT(vm.count("command"));

        hpx::id_type target = hpx::naming::get_id_from_locality_id
            (vm["target"].as<boost::uint32_t>());

        std::string raw_cmd = vm["command"].as<std::string>();

        command_interpreter(target, raw_cmd);
    }

    hpx::disconnect();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    using boost::program_options::value;
    using boost::program_options::options_description;

    // Configure application-specific options.
    options_description cmdline("Usage: " HPX_APPLICATION_STRING " [options]");

    cmdline.add_options()
        ( "target", value<boost::uint32_t>()->default_value(0)
        , "locality to connect to")
        ( "command,c", value<std::string>(), "command to execute")
    ;

    // Disable loading of all external components.
    std::vector<std::string> cfg;
    cfg.push_back("hpx.components.load_external=0");
    HPX_STD_FUNCTION<void()> empty;

    return hpx::init(cmdline, argc, argv, cfg, empty, empty,
        hpx::runtime_mode_connect);
}

