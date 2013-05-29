////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_ED5406AF_2905_4BDD_B8A4_78D28B9E38AC)
#define OCTOPUS_ED5406AF_2905_4BDD_B8A4_78D28B9E38AC

#include <octopus/octree/octree_server.hpp>

#include <boost/cstdint.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/version.hpp>

#define OCTOPUS_WRITER_BASE_VERSION 0x01

namespace octopus
{

struct OCTOPUS_EXPORT writer_base
{
    virtual ~writer_base() {}

    // FIXME: Would be nice to get rid of the "time" parameter.
    virtual void begin_epoch(
        octree_server& e
      , double time
      , std::string const& file
        )
    { }

    virtual void end_epoch(octree_server& e) { }

    virtual void operator()(octree_server& e) = 0;

    virtual writer_base* clone() const = 0;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {}
};

struct writer
{
  private:
    boost::shared_ptr<writer_base> ptr_;

  public:
    writer() : ptr_() {}

    writer(writer const& other)
      : ptr_(other.ptr_)
    {}

    writer& operator=(writer const& other)
    {
        ptr_ = other.ptr_;
        return *this;
    }

    writer& operator=(writer_base const& other)
    {
        ptr_.reset(other.clone());
        return *this;
    }

    void begin_epoch(
        octree_server& e
      , double time
      , std::string const& file
        )
    {
        ptr_->begin_epoch(e, time, file);
    }

    void end_epoch(octree_server& e) const
    {
        ptr_->end_epoch(e);
    }

    template <typename Derived>
    Derived* cast() const
    {
        return dynamic_cast<Derived*>(ptr_.get());
    }

    void operator()(octree_server& e) const
    {
        (*ptr_)(e);
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        bool has_ptr = ptr_;

        ar & has_ptr;

        if (has_ptr)
        {
            writer_base const* r = ptr_.get();
            ar & r;
        } 
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        bool has_ptr = false;

        ar & has_ptr;

        if (has_ptr)
        {
            writer_base* r = 0;
            ar & r;
            ptr_.reset(r);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

}

BOOST_CLASS_VERSION(octopus::writer_base, OCTOPUS_WRITER_BASE_VERSION)
BOOST_CLASS_TRACKING(octopus::writer_base, boost::serialization::track_never)


#endif // OCTOPUS_ED5406AF_2905_4BDD_B8A4_78D28B9E38AC

