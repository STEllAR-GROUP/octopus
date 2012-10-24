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
    boost::uint64_t variable_index_;
    std::string variable_name_;
    std::string file_name_;
    bool merged_;

    void open_locked(boost::uint64_t step, mutex_type::scoped_lock& l);

    void close_locked(mutex_type::scoped_lock& l);

    void merge_locked(mutex_type::scoped_lock& l);

  public:
    single_variable_silo_writer()
      : mtx_()
      , file_()
      , directory_names_()
      , step_(0)
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
      , std::string const& file_name = "U_S%06u_L%06u.silo"
        )
      : mtx_()
      , file_(0)
      , directory_names_()
      , step_(0)
      , variable_index_(variable_index)
      , variable_name_(variable_name)
      , file_name_(file_name)
      , merged_(false)
    {}

    ~single_variable_silo_writer()
    {
        mutex_type::scoped_lock l(mtx_);
        merge_locked(l);
        close_locked(l);
    }

    void open(boost::uint64_t step)
    {
        mutex_type::scoped_lock l(mtx_);
        open_locked(step, l);
    }

    void close()
    {
        mutex_type::scoped_lock l(mtx_);
        close_locked(l);
    }

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

