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
    result_type operator()(std::string const& key, A& data, B dflt) const
    {
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

// }}}
 
///////////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(
    std::ostream& os
  , config_data const& cfg
    )
{
    char const* fmt = "%1% %|40t| = %2%";

    // NOTE: Last item should not have a newline after it.
    os
        << (boost::format(fmt) % "octopus.dimensional_size"
                               % cfg.dimensional_size) << "\n"
        << (boost::format(fmt) % "octopus.temporal_prediction_gap"
                               % cfg.temporal_prediction_gap)
    ;

    return os;
}

config_data config_from_ini()
{
    config_data cfg;

    config_reader reader;

    reader
        ("octopus.dimensional_size", cfg.dimensional_size, 12) 
        ("octopus.temporal_prediction_gap", cfg.temporal_prediction_gap, 10) 
    ;

    return cfg;
}

}

