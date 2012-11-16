////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_3BF02C42_6BDB_4469_B6D1_8B8F393C63A5)
#define OCTOPUS_3BF02C42_6BDB_4469_B6D1_8B8F393C63A5

#include <octopus/io/writer.hpp>

#include <hpx/util/base_object.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/version.hpp>

#include <string>

#include <silo.h>

// FIXME: If copied with an open file, this should probably assert.
// FIXME: Naming for silo_writer breaks the general convention.

namespace octopus
{ 

struct OCTOPUS_EXPORT single_variable_silo_writer : writer_base
{
  private:
    typedef hpx::lcos::local::mutex mutex_type;

    mutex_type mtx_;
    DBfile* file_;

    std::vector<std::string> directory_names_;
    boost::uint64_t step_;
    double time_;
    boost::uint64_t variable_index_;
    std::string variable_name_;
    std::string file_name_;
    bool merged_;

    void start_write_locked(
        boost::uint64_t step
      , double time
      , std::string const& file 
      , mutex_type::scoped_lock& l
        );

    void stop_write_locked(mutex_type::scoped_lock& l);

    void merge_locked(mutex_type::scoped_lock& l);

  public:
    single_variable_silo_writer()
      : mtx_()
      , file_()
      , directory_names_()
      , step_(0)
      , time_(0.0)
      , variable_index_()
      , variable_name_()
      , file_name_()
      , merged_(false)
    {}

    single_variable_silo_writer(
        single_variable_silo_writer const& other
        )
      : mtx_()
      , file_()
      , directory_names_(other.directory_names_)
      , step_(other.step_)
      , time_(other.time_)
      , variable_index_(other.variable_index_)
      , variable_name_(other.variable_name_)
      , file_name_(other.file_name_)
      , merged_(false)
    {}
 
    single_variable_silo_writer(
        boost::uint64_t variable_index
      , std::string const& variable_name 
                                     // These should be sufficient default
                                     // widths (for now).
      , std::string const& file_name = "U_L%06u_S%06u.silo"
        )
      : mtx_()
      , file_(0)
      , directory_names_()
      , step_(0)
      , time_(0.0)
      , variable_index_(variable_index)
      , variable_name_(variable_name)
      , file_name_(file_name)
      , merged_(false)
    {}

    ~single_variable_silo_writer()
    {
        mutex_type::scoped_lock l(mtx_);
        merge_locked(l);
        stop_write_locked(l);
    }

    void start_write(boost::uint64_t step, double time, std::string const& file)
    {
        mutex_type::scoped_lock l(mtx_);
        start_write_locked(step, time, file, l);
    }

    void start_write(boost::uint64_t step, double time)
    {
        mutex_type::scoped_lock l(mtx_);
        start_write_locked(step, time, file_name_, l);
    }

    void stop_write()
    {
        mutex_type::scoped_lock l(mtx_);
        stop_write_locked(l);
    }

    void begin_epoch(octree_server& e, std::string const& file);

    void begin_epoch(octree_server& e)
    {
        begin_epoch(e, file_name_);
    }

    void end_epoch(octree_server& e);

    void merge()
    {
        mutex_type::scoped_lock l(mtx_);
        merge_locked(l);
    }

    // IMPLEMENT: Do actual I/O in a separate OS-thread.
    void operator()(octree_server& e);

    writer_base* clone() const
    {
        return new single_variable_silo_writer(variable_index_
                                             , variable_name_
                                             , file_name_);
    } 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & hpx::util::base_object_nonvirt<writer_base>(*this);
        ar & step_;
        ar & time_;
        ar & variable_index_;
        ar & variable_name_;
        ar & file_name_;
    }
};

}

BOOST_CLASS_EXPORT_GUID(octopus::single_variable_silo_writer
                      , "single_variable_silo_writer")
BOOST_CLASS_TRACKING(octopus::single_variable_silo_writer
                   , boost::serialization::track_never)

#endif // OCTOPUS_3BF02C42_6BDB_4469_B6D1_8B8F393C63A5

