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

// FIXME FIXME FIXME
// The entire command string has to be parsed with cgi_query_to_ini, including
// the command. We can't do a pop_front on std::vector, so we need to handle
// the first cgi query (which should be the command) separately from the rest. 
//
// Also, update_config should take any number of arguments greater than 1.
// FIXME FIXME FIXME

typedef std::pair<std::string, std::vector<std::string> > command;

command parse_cgi_query_to_command(std::string raw_cmd)
{
    std::size_t found = raw_cmd.find_first_of("&");

    if (found == std::string::npos)
        return std::pair<std::string, std::vector<std::string> >
            (raw_cmd, std::vector<std::string>());

    std::pair<std::string, std::vector<std::string> > parsed_cmd;

    parsed_cmd.first = raw_cmd.substr(0, found);
    raw_cmd = raw_cmd.substr(found + 1);

    boost::algorithm::split(parsed_cmd.second, raw_cmd,
        boost::algorithm::is_any_of("&"),
        boost::algorithm::token_compress_on);

    return parsed_cmd;
}

///////////////////////////////////////////////////////////////////////////////
void command_interpreter(hpx::id_type const& target, std::string const& raw_cmd)
{
    command parsed_cmd = parse_cgi_query_to_command(raw_cmd);

    std::string& cmd = parsed_cmd.first;
    std::vector<std::string>& args = parsed_cmd.second;

    if (cmd.empty())
    {
        HPX_THROW_EXCEPTION(hpx::bad_parameter, "command_interpreter",
                "no command specified");
        return;
    }

    // jQuery adds a timestep to HTTP GET requests in the form of _= as a cache
    // buster.

    if (  args.back().size() >= 2
       && (args.back()[0] == '_' && args.back()[1] == '='))
        args.resize(args.size() - 1);

    ///////////////////////////////////////////////////////////////////////////
    // Update the INI config. 
    if (cmd == "update_config")
    {
        if (args.size() == 0)
        {
            HPX_THROW_EXCEPTION(hpx::bad_parameter, "command_interpreter",
                "get_config takes 1 or more arguments");
            return;
        }

        hpx::util::section ini;
        ini.parse("runtime update", args, false);

        hpx::async<update_live_config_action>(target, ini).get();

        std::cout << ini_to_json(ini) << std::flush;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Gets the current config
    else if (cmd == "get_config")
    {
        if (args.size() != 1)
        {
            HPX_THROW_EXCEPTION(hpx::bad_parameter, "command_interpreter",
                "get_config takes only 1 argument");
            return;
        }

        hpx::util::section ini
            = hpx::async<get_live_config_action>(target, args[0]).get();

        std::cout << ini_to_json(ini) << std::flush;
    }

    else
        {
            HPX_THROW_EXCEPTION(hpx::bad_parameter, "command_interpreter",
                "unknown command");
            return;
        }
}

///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map& vm)
{
    try
    {
        if (!vm.count("command"))
        {
            HPX_THROW_EXCEPTION(hpx::bad_parameter, "hpx_main",
                "no command specified");
            return 1;
        }

        hpx::id_type target = hpx::naming::get_id_from_locality_id
            (vm["target"].as<boost::uint32_t>());

        std::string raw_cmd = vm["command"].as<std::string>();

        command_interpreter(target, raw_cmd);
    }

    catch (hpx::exception const& e)
    {
        std::cout << "{what}: "        << hpx::get_error_what(e) << "\n";
        std::cout << "{function}: "    << hpx::get_error_function_name(e) << "\n";
        std::cout << "{file}: "        << hpx::get_error_file_name(e) << "\n";
        std::cout << "{line}: "        << hpx::get_error_line_number(e) << "\n";
        std::cout << "{locality-id}: " << hpx::get_error_locality_id(e) << "\n";
        std::cout << "{os-thread}: "   << hpx::get_error_os_thread(e) << "\n";
        std::cout << "{thread-id}: "   << std::hex << hpx::get_error_thread_id(e) << "\n";
        //std::cout << "{stack-trace}: " << hpx::get_error_backtrace(e) << "\n";
        std::cout << std::flush;

        hpx::disconnect();
        return 1;
    }  

    catch (boost::system::system_error const& e) {
        std::cout << "{what}: " << e.what() << "\n";
        std::cout << std::flush;        
        hpx::disconnect();
        return 1;
    }
    catch (std::exception const& e) {
        std::cout << "{what}: " << e.what() << "\n";
        std::cout << std::flush;
        hpx::disconnect();
        return 1;
    }
    catch (...) {
        std::cout << "{what}: unknown exception\n";
        std::cout << std::flush;
        hpx::disconnect();
        return 1;
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
        ( "command", value<std::string>(), "command to execute")
    ;

    // Disable loading of all external components.
    std::vector<std::string> cfg;
    cfg.push_back("hpx.components.load_external=0");
    HPX_STD_FUNCTION<void()> empty;

    return hpx::init(cmdline, argc, argv, cfg, empty, empty,
        hpx::runtime_mode_connect);
}

