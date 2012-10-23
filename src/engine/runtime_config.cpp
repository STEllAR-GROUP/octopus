////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_fwd.hpp>
#include <hpx/exception.hpp>
#include <hpx/runtime.hpp>

#include <octopus/engine/runtime_config.hpp>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <iomanip>


namespace octopus
{

///////////////////////////////////////////////////////////////////////////////
// {{{ Configuration utilities

inline bool has_config_entry(std::string const& key)
{                                                                            
    if (NULL == hpx::get_runtime_ptr())                                           
        return false;                                                           
    return hpx::get_runtime().get_config().has_entry(key);                  
}  

struct config_reader
{
    typedef config_reader const& result_type;

    template <typename A, typename B>
    result_type operator()(std::string const& param, A& data, B dflt) const
    {
        std::string key("octopus.");
        key += param;

        try
        {
            if (has_config_entry(key))
            {
                data = boost::lexical_cast<A>(hpx::get_config_entry(key, ""));
            }
            else
            {
                data = dflt;
            }
    
        }
    
        catch (boost::bad_lexical_cast&)
        {
            std::string msg = boost::str(boost::format(
                "bad INI parameter, '%1%' is not a valid value for %2%")
                % hpx::get_config_entry(key, "") % key);
            HPX_THROW_EXCEPTION(hpx::bad_parameter,
                "octopus::config_from_ini", msg);
        }

        return *this;
    }
};

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
    #define OCTOPUS_FORMAT_OPTION(option)                                   \
        format_option("octopus." BOOST_PP_STRINGIZE(option), cfg.option)    \
        /**/

    // NOTE: Last item should not have a newline after it.
    os
        << OCTOPUS_FORMAT_OPTION(max_refinement_level) << "\n"
        << OCTOPUS_FORMAT_OPTION(grid_node_length) << "\n"
        << OCTOPUS_FORMAT_OPTION(spatial_domain) << "\n"
        << OCTOPUS_FORMAT_OPTION(initial_spatial_step) << "\n"
        << OCTOPUS_FORMAT_OPTION(runge_kutta_order) << "\n"
        << OCTOPUS_FORMAT_OPTION(x_reflect) << "\n"
        << OCTOPUS_FORMAT_OPTION(y_reflect) << "\n"
        << OCTOPUS_FORMAT_OPTION(z_reflect) << "\n"
        << OCTOPUS_FORMAT_OPTION(temporal_prediction_gap)
    ;

    #undef OCTOPUS_FORMAT_OPTION

    return os;
}

config_data config_from_ini()
{
    config_data cfg;

    config_reader reader;

    // FIXME: Math in INI would make this smoother, some of these settings
    // should default to a formula not a hard-coded value.
    reader
        ("max_refinement_level", cfg.max_refinement_level, 0) 
        ("grid_node_length", cfg.grid_node_length, 12) 
        ("spatial_domain", cfg.spatial_domain, 1.5e-4) 
        ("initial_spatial_step", cfg.initial_spatial_step, 3.75e-5) 
        ("runge_kutta_order", cfg.runge_kutta_order, 1) 
        ("x_reflect", cfg.x_reflect, false) 
        ("y_reflect", cfg.y_reflect, false) 
        ("z_reflect", cfg.z_reflect, true) 
        ("temporal_prediction_gap", cfg.temporal_prediction_gap, 10) 
    ;

    return cfg;
}

}

