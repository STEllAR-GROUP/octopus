////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/visit/visit_simulation_server.hpp>

#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>

namespace octopus
{

void visit_simulation_server::start(
    std::string const& name
  , std::string const& sim_file
  , std::string const& exe
  , std::vector<std::string> const& args
    // environment is a std::map<std::string, std::string>
  , boost::process::environment env // by value for swap
    )
{
    namespace bp = boost::process;

    OCTOPUS_ASSERT(!sim_);

    VisItSetupEnvironment();

    VisItInitializeSocketAndDumpSimFile(
        name.c_str() // simulation name
      , NULL // description
      , NULL // simulation directory, leaving it empty appears to be fine
      , NULL // unknown, called inputfile in visit headers
      , NULL // unknown, called guifile in visit headers
      , sim_file.c_str() // simulation file
        ); 

    bp::context ctx;
    ctx.environment = bp::self::get_environment(); 

    // Merge the current environment and the requested environment (complexity
    // ugh), overwriting variables in the current environment with variables
    // specified in the env argument.
    env.insert(ctx.environment.begin(), ctx.environment.end());
    env.swap(ctx.environment);

    sim_.reset(new bp::child(bp::launch(exe, args, ctx)));

    VisItDetectInput(0, -1);
     
    if (!VisItAttemptToCompleteConnection()) 
        OCTOPUS_ASSERT_MSG(false, "visit did not connect");
}

void visit_simulation_server::evaluate(std::string const& source)
{
    OCTOPUS_ASSERT(sim_);

    VisItExecuteCommand(source.c_str());    
}

// FIXME: I don't really do what I'm suppose to :(.
void visit_simulation_server::terminate()
{
    OCTOPUS_ASSERT(sim_);

    VisItExecuteCommand("Close()");
    VisItExecuteCommand("quit()");
}

}

