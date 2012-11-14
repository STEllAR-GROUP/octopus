////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_3D61DA0E_286D_486B_A4D5_BE452310CC5D)
#define OCTOPUS_3D61DA0E_286D_486B_A4D5_BE452310CC5D

#include <octopus/visit/visit_simulation_server.hpp>
#include <octopus/filesystem.hpp>

#include <hpx/runtime/components/client_base.hpp>
#include <hpx/runtime/components/stubs/stub_base.hpp>

namespace octopus
{

struct visit_simulation_client
  : hpx::components::client_base<
        visit_simulation_client
      , hpx::components::stub_base<visit_simulation_server>
    >
{
    typedef hpx::components::client_base<
        visit_simulation_client
      , hpx::components::stub_base<visit_simulation_server>
    > base_type;

    visit_simulation_client() : base_type() {}

    visit_simulation_client(visit_simulation_client const& other)
      : base_type(other.gid_)
    {}

    visit_simulation_client(BOOST_RV_REF(visit_simulation_client) other)
      : base_type(boost::move(other.gid_))
    {}

    visit_simulation_client& operator=(
        BOOST_COPY_ASSIGN_REF(visit_simulation_client) other
        )
    {
        this->gid_ = other.gid_;
        return *this;
    }

    visit_simulation_client& operator=(
        BOOST_RV_REF(visit_simulation_client) other
        )
    {
        this->gid_ = boost::move(other.gid_);
        return *this;
    }

    ///////////////////////////////////////////////////////////////////////////
    // {{{ start
    hpx::future<void> start_async(
        std::string const& name
        )
    {
        std::string sim_file = join_paths(current_path(), name + ".sim2");
        return start_async(name, sim_file);
    }

    hpx::future<void> start_async(
        std::string const& name
      , std::string const& sim_file
        )
    {
        std::string exe = join_paths(OCTOPUS_VISIT_DRIVER_ROOT, "bin", "visit"); 
        return start_async(name, sim_file, exe);
    }

    hpx::future<void> start_async(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exe
        )
    {
//        std::vector<std::string> args { "", "-fullscreen", "-cli", "-o", sim_file };
        std::vector<std::string> args { "", "-o", sim_file };
        return start_async(name, sim_file, exe, args);  
    }

    hpx::future<void> start_async(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exe
      , std::vector<std::string> const& args
        )
    {
        boost::process::environment env;
        return start_async(name, sim_file, exe, args, env);  
    }

    hpx::future<void> start_async(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exe
      , std::vector<std::string> const& args
      , boost::process::environment const& env 
        )
    {
        return hpx::async<visit_simulation_server::start_action>
            (get_gid(), name, sim_file, exe, args, env);
    }

    ///////////////////////////////////////////////////////////////////////////
    void start(
        std::string const& name
        )
    {
        start_async(name).get();
    }

    void start(
        std::string const& name
      , std::string const& sim_file
        )
    {
        start_async(name, sim_file).get();
    }

    void start(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exe
        )
    {
        start_async(name, sim_file, exe).get();
    }

    void start(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exe
      , std::vector<std::string> const& args
        )
    {
        start_async(name, sim_file, exe, args).get();
    }

    void start(
        std::string const& name
      , std::string const& sim_file
      , std::string const& exe
      , std::vector<std::string> const& args
      , boost::process::environment const& env 
        )
    {
        start_async(name, sim_file, exe, args, env).get();
    }
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ terminate
    hpx::future<void> terminate_async()
    {
        return hpx::async<visit_simulation_server::terminate_action>
            (get_gid());
    } 

    void terminate()
    {
        terminate_async().get();
    } 
    // }}}

    ///////////////////////////////////////////////////////////////////////////
    // {{{ evaluate
    hpx::future<void> evaluate_async(std::string const& source)
    {
        return hpx::async<visit_simulation_server::evaluate_action>
            (get_gid(), source);
    } 

    void evaluate(std::string const& source)
    {
        evaluate_async(source).get();
    } 

    void operator<<(std::string const& source)
    {
        evaluate_async(source).get();
    } 
    // }}}
};

}

#endif // OCTOPUS_3D61DA0E_286D_486B_A4D5_BE452310CC5D

