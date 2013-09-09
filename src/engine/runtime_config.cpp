////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/engine/ini.hpp>
#include <octopus/engine/runtime_config.hpp>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <iomanip>

namespace octopus
{

// TODO: Refine this stuff and move it into a header.

///////////////////////////////////////////////////////////////////////////////
// {{{ Formatting utilities
boost::format format_option(
    std::string const& name
  , bool value
    )
{
    char const* const fmt = "%1% %|40t| = %2%";
    return ( boost::format(fmt) % name
           % boost::io::group(std::boolalpha, value));
}

boost::format format_option(
    std::string const& name
  , double value
    )
{
    char const* const fmt = "%1% %|40t| = %2%";
    return ( boost::format(fmt) % name
           % boost::io::group(std::scientific, value));
}

template <typename T>
boost::format format_option(
    std::string const& name
  , T const& value
    )
{
    char const* const fmt = "%1% %|40t| = %2%";
    return (boost::format(fmt) % name % value);
}

// }}}
 
///////////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(
    std::ostream& os
  , config_data const& cfg
    )
{
    #define OCTOPUS_FORMAT_OPTION(option)                        \
        format_option(BOOST_PP_STRINGIZE(option), cfg.option)    \
        /**/

    // NOTE: Last item should not have a newline after it.
    os
        << "[octopus]\n"
        << OCTOPUS_FORMAT_OPTION(levels_of_refinement) << "\n"

        << OCTOPUS_FORMAT_OPTION(runge_kutta_order) << "\n"
        << OCTOPUS_FORMAT_OPTION(reflect_on_z) << "\n"

        << OCTOPUS_FORMAT_OPTION(spatial_domain) << "\n"
        << OCTOPUS_FORMAT_OPTION(grid_node_length) << "\n"

        << OCTOPUS_FORMAT_OPTION(temporal_domain) << "\n"
        << OCTOPUS_FORMAT_OPTION(temporal_prediction_gap) << "\n"

        << OCTOPUS_FORMAT_OPTION(output_frequency) << "\n"

        << OCTOPUS_FORMAT_OPTION(checkpoint_file) << "\n"
        << OCTOPUS_FORMAT_OPTION(load_checkpoint)
    ;

    #undef OCTOPUS_FORMAT_OPTION

    return os;
}

config_data config_from_ini()
{
    config_data cfg;

    config_reader reader("octopus");

    // FIXME: Math in INI would make this smoother, some of these settings
    // should default to a formula not a hard-coded value.
    reader
        ("levels_of_refinement", cfg.levels_of_refinement, 0) 

        ("runge_kutta_order", cfg.runge_kutta_order, 3) 
        ("reflect_on_z", cfg.reflect_on_z, true) 

        ("spatial_domain", cfg.spatial_domain, 1.5) 
        ("grid_node_length", cfg.grid_node_length, 14) // FIXME: This should be
                                                       // computed dynamically. 

        ("temporal_domain", cfg.temporal_domain, 1.0) 
        ("temporal_prediction_gap", cfg.temporal_prediction_gap, 10) 
        
        ("output_frequency", cfg.output_frequency, 0.005)

        ("checkpoint_file", cfg.checkpoint_file, "checkpoint_L%06u.bin")
        ("load_checkpoint", cfg.load_checkpoint, false)
    ;

    return cfg;
}

}

