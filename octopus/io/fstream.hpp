////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_A773E6B8_8BAB_4336_8503_87049CFA0C1B)
#define OCTOPUS_A773E6B8_8BAB_4336_8503_87049CFA0C1B

#include <octopus/io/writer.hpp>

#include <hpx/util/base_object.hpp>

#include <fstream>

// FIXME: If copied with an open file, this should probably assert.
// FIXME: The mutex is bad here, we end up serializing all I/O and leaving a ton
// of threads waiting on the stack.

namespace octopus
{ 

typedef hpx::util::function<
    void(
        octree_server&
      , std::ofstream&
        )
> fstream_writer_function;

struct OCTOPUS_EXPORT fstream_writer : writer_base
{
  private:
    typedef hpx::lcos::local::mutex mutex_type;

    mutex_type mtx_;
    std::ofstream file_;
    std::string file_name_;
    fstream_writer_function f_;
    double time_;

  public:
    fstream_writer()
      : mtx_()
      , file_()
      , file_name_()
      , f_()
      , time_(0.0)
    {}

    fstream_writer(
        fstream_writer const& other
        )
      : mtx_()
      , file_()
      , file_name_(other.file_name_)
      , f_(other.f_)
      , time_(other.time_)
    {}
 
    fstream_writer(
        fstream_writer_function const& f 
      , std::string const& file_name
        )
      : mtx_()
      , file_()
      , file_name_(file_name)
      , f_(f)
      , time_()
    {}

    ~fstream_writer()
    {
        stop_write();
    }

    void start_write(
        boost::uint64_t step
      , double time
      , std::string const& file
        );

    void stop_write();

    void begin_epoch(octree_server& e, double time, std::string const& file);

    void end_epoch(octree_server& e);

    // IMPLEMENT: Do actual I/O in a separate OS-thread.
    void operator()(octree_server& e);

    writer_base* clone() const
    {
        return new fstream_writer(f_, file_name_);
    } 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & hpx::util::base_object_nonvirt<writer_base>(*this);
        ar & file_name_;
        ar & f_;
        ar & time_;
    }
};

}

BOOST_CLASS_EXPORT_GUID(octopus::fstream_writer, "fstream_writer")
BOOST_CLASS_TRACKING(octopus::fstream_writer, boost::serialization::track_never)

#endif // OCTOPUS_A773E6B8_8BAB_4336_8503_87049CFA0C1B

