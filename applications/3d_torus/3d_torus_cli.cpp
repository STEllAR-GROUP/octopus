////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_init.hpp>
#include <hpx/include/plain_actions.hpp>
#include <hpx/async.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

// Dummy.
void set_zoom(std::string const& arg) { BOOST_ASSERT(false); }
HPX_PLAIN_ACTION(set_zoom, set_zoom_action);
HPX_ACTION_HAS_CRITICAL_PRIORITY(set_zoom_action);

// Dummy.
void set_angle(std::string const& arg1, std::string const& arg2)
{ BOOST_ASSERT(false); }
HPX_PLAIN_ACTION(set_angle, set_angle_action);
HPX_ACTION_HAS_CRITICAL_PRIORITY(set_angle_action);

// Dummy.
void update_kappa(double k) { BOOST_ASSERT(false); }
HPX_PLAIN_ACTION(update_kappa, update_kappa_action);
HPX_ACTION_HAS_CRITICAL_PRIORITY(update_kappa_action);

char const* const help = "commands: set_zoom [double], "
                                   "set_angle [int, int], "
                                   "set_kappa [double]";

///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map& vm)
{
    {
        BOOST_ASSERT(vm.count("target"));

        hpx::id_type target = hpx::naming::get_id_from_locality_id
            (vm["target"].as<boost::uint32_t>());

        // Print out the available commands.
        std::cout << help << std::endl << "> ";

        // Enter the interpreter loop.
        std::string line;
        while (std::getline(std::cin, line))
        {
            boost::algorithm::trim(line);

            std::vector<std::string> cmd;
            boost::algorithm::split(cmd, line,
                boost::algorithm::is_any_of(" \t\n"),
                boost::algorithm::token_compress_on);

            if (!cmd.empty() && !cmd[0].empty()) 
            {
                if (cmd[0] == "set_zoom")
                {
                    if (cmd.size() == 2)
                    {
                        hpx::async<set_zoom_action>(target, cmd[1]).get();
                    }
                    else
                    {
                        std::cout << "error: invalid command '"
                                  << line << "'" << std::endl
                                  << help << std::endl;
                    }
                }
                else if (cmd[0] == "set_angle")
                {
                    if (cmd.size() == 3)
                    {
                        hpx::async<set_angle_action>
                            (target, cmd[1], cmd[2]).get();
                    }
                    else
                    {
                        std::cout << "error: invalid command '"
                                  << line << "'" << std::endl
                                  << help << std::endl;
                    }
                }
                else if (cmd[0] == "set_kappa")
                {
                    if (cmd.size() == 2)
                    {
                        hpx::async<update_kappa_action>(target
                            , boost::lexical_cast<double>(cmd[1])).get();
                    }
                    else
                    {
                        std::cout << "error: invalid command '"
                                  << line << "'" << std::endl
                                  << help << std::endl;
                    }
                }
                else
                {
                    std::cout << "error: invalid command '"
                              << line << "'" << std::endl
                              << help << std::endl;
                }
            }

            std:: cout << "> ";
        }
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
        ("target", value<boost::uint32_t>(), "locality to connect to")
    ;

    // Disable loading of all external components.
    std::vector<std::string> cfg;
    cfg.push_back("hpx.components.load_external=0");
    HPX_STD_FUNCTION<void()> empty;

    return hpx::init(cmdline, argc, argv, cfg, empty, empty,
        hpx::runtime_mode_connect);
}

