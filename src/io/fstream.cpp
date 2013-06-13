////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/include/actions.hpp>
#include <hpx/lcos/future_wait.hpp>

#include <octopus/io/fstream.hpp>
#include <octopus/engine/engine_interface.hpp>

namespace octopus
{

void fstream_writer::start_write(
    boost::uint64_t step
  , double time
    )
{
    mutex_type::scoped_lock l(mtx_);

    // Make sure we closed the last epoch.
    if (file_.is_open())
        file_.close();

    try
    {
        std::string s = boost::str( boost::format(file_name_)
                                  % hpx::get_locality_id() % step);
        file_.open(s, std::fstream::binary); 
    }
    // FIXME: Catch the specific boost.format exception.
    catch (...)
    {
        try
        {
            std::string s = boost::str( boost::format(file_name_)
                                      % hpx::get_locality_id());
            file_.open(s, std::fstream::binary); 
        }
        // FIXME: Catch the specific boost.format exception.
        catch (...)
        {
            file_.open(file_name_, std::fstream::binary); 
        }
    }

    OCTOPUS_ASSERT(file_.is_open());
}

void fstream_writer::stop_write()
{
    mutex_type::scoped_lock l(mtx_);

    if (file_.is_open())
        file_.close();
}

void fstream_perform_start_write(
    boost::uint64_t step
  , double time
    )
{
    science().output.cast<fstream_writer>()->start_write(step, time);
}

void fstream_perform_stop_write()
{
    science().output.cast<fstream_writer>()->stop_write();
}

}

HPX_PLAIN_ACTION(octopus::fstream_perform_start_write
               , fstream_perform_start_write_action);
HPX_PLAIN_ACTION(octopus::fstream_perform_stop_write
               , fstream_perform_stop_write_action);

namespace octopus
{

void fstream_writer::begin_epoch(
    octree_server& e
  , double time
    )
{
    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    fstream_perform_start_write_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
    {
        futures.emplace_back(hpx::async(act, targets[i], e.get_step(), time));
    }

    hpx::wait(futures);
}

void fstream_writer::end_epoch(octree_server& e)
{
    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    fstream_perform_stop_write_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
        futures.push_back(hpx::async(act, targets[i]));

    hpx::wait(futures);
}

void fstream_writer::operator()(octree_server& e)
{
    mutex_type::scoped_lock l(mtx_);

    OCTOPUS_ASSERT(f_);
    OCTOPUS_ASSERT(file_.is_open());

    f_(e, file_);
}

}

