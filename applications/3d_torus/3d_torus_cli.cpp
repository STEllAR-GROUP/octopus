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
void set_momentum_conservation(std::string const& arg) { BOOST_ASSERT(false); }
HPX_PLAIN_ACTION(set_momentum_conservation, set_momentum_conservation_action);

///////////////////////////////////////////////////////////////////////////////
int hpx_main(boost::program_options::variables_map& vm)
{
    {
        BOOST_ASSERT(vm.count("target"));
        BOOST_ASSERT(vm.count("set-momentum-conservation"));

        hpx::id_type target = hpx::naming::get_id_from_locality_id
            (vm["target"].as<boost::uint32_t>());

        std::string scheme = vm["set-momentum-conservation"].as<std::string>();
        hpx::async<set_momentum_conservation_action>(target, scheme).get();
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
        ("target"
        , value<boost::uint32_t>()->default_value(0)
        , "locality to connect to")
        ( "set-momentum-conservation"
        , value<std::string>()
        , "new momentum conservation scheme")
    ;

    // Disable loading of all external components.
    std::vector<std::string> cfg;
    cfg.push_back("hpx.components.load_external=0");
    HPX_STD_FUNCTION<void()> empty;

    return hpx::init(cmdline, argc, argv, cfg, empty, empty,
        hpx::runtime_mode_connect);
}

