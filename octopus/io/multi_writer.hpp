////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_F255A1F0_667F_4FED_8AD0_495BF1DD0284)
#define OCTOPUS_F255A1F0_667F_4FED_8AD0_495BF1DD0284

#include <octopus/io/writer.hpp>

#include <boost/serialization/vector.hpp>

namespace octopus
{

struct multi_writer : writer_base
{
  private:
    std::vector<writer> writers_;

  public:
    multi_writer() : writers_() {}

    multi_writer(multi_writer const& other)
      : writers_(other.writers_)
    {}

    multi_writer(BOOST_RV_REF(multi_writer) other)
      : writers_(other.writers_)
    {}

    multi_writer(std::vector<writer> const& writers)
      : writers_(writers)
    {}

    multi_writer(BOOST_RV_REF(std::vector<writer>) writers)
      : writers_(writers)
    {}

    void add_writer(writer const& w)
    {
        writers_.push_back(w);  
    }

    void add_writer(BOOST_RV_REF(writer) w)
    {
        writers_.push_back(w);  
    }

    void begin_epoch(octree_server& e, double time)
    {
        for (boost::uint64_t i = 0; i < writers_.size(); ++i)
            writers_[i].begin_epoch(e, time);
    }

    void end_epoch(octree_server& e)
    {
        for (boost::uint64_t i = 0; i < writers_.size(); ++i)
            writers_[i].end_epoch(e);
    }

    // IMPLEMENT: Do actual I/O in a separate OS-thread.
    void operator()(octree_server& e)
    {
        for (boost::uint64_t i = 0; i < writers_.size(); ++i)
            writers_[i](e);
    }

    writer_base* clone() const
    {
        return new multi_writer(writers_);
    } 

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & hpx::util::base_object_nonvirt<writer_base>(*this);
        ar & writers_;
    }
};

}

BOOST_CLASS_EXPORT_GUID(octopus::multi_writer, "multi_writer")
BOOST_CLASS_TRACKING(octopus::multi_writer, boost::serialization::track_never)

#endif // OCTOPUS_F255A1F0_667F_4FED_8AD0_495BF1DD0284

