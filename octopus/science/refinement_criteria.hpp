////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_95F7B0D1_2246_4B89_9B80_BAC4C4F22C0E)
#define OCTOPUS_95F7B0D1_2246_4B89_9B80_BAC4C4F22C0E

#include <octopus/octree/octree_server.hpp>

#include <boost/cstdint.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace octopus
{

struct OCTOPUS_EXPORT refinement_criteria_base
{
    virtual ~refinement_criteria_base() {}

    virtual bool refine(octree_server& e, child_index idx) = 0;

    virtual bool unrefine(octree_server& e, child_index idx) = 0; 

    virtual refinement_criteria_base* clone() const = 0;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {}
};

struct refinement_criteria
{
  private:
    boost::shared_ptr<refinement_criteria_base> ptr_;

  public:
    refinement_criteria() : ptr_() {}

    refinement_criteria(refinement_criteria const& other)
      : ptr_(other.ptr_)
    {}

    refinement_criteria& operator=(refinement_criteria const& other)
    {
        ptr_ = other.ptr_;
        return *this;
    }

    refinement_criteria& operator=(refinement_criteria_base const& other)
    {
        ptr_.reset(other.clone());
        return *this;
    }

    bool refine(octree_server& e, child_index idx) const
    {
        return ptr_->refine(e, idx);
    }

    bool unrefine(octree_server& e, child_index idx) const 
    {
        return ptr_->unrefine(e, idx);
    }

    template <typename Derived>
    Derived* cast() const
    {
        return dynamic_cast<Derived*>(ptr_.get());
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        bool has_ptr = ptr_;

        ar & has_ptr;

        if (has_ptr)
        {
            refinement_criteria_base const* r = ptr_.get();
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
            refinement_criteria_base* r = 0;
            ar & r;
            ptr_.reset(r);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

}

#endif // OCTOPUS_95F7B0D1_2246_4B89_9B80_BAC4C4F22C0E

