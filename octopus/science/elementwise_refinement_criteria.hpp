////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_EC4E7955_BF33_473F_9A44_FF274B50EEE8)
#define OCTOPUS_EC4E7955_BF33_473F_9A44_FF274B50EEE8

#include <hpx/util/base_object.hpp>
#include <hpx/util/detail/serialization_registration.hpp>

#include <octopus/engine/engine_interface.hpp>
#include <octopus/science/refinement_criteria.hpp>

namespace octopus
{

template <typename Derived>
struct OCTOPUS_EXPORT elementwise_refinement_criteria_base
  : refinement_criteria_base
{
    Derived& derived()
    {
        return static_cast<Derived&>(*this);
    }

    bool refine(octree_server& e, child_index idx)
    {
        boost::uint64_t const bw = science().ghost_zone_width;
        boost::uint64_t const gnx = config().grid_node_length;
    
        for (boost::uint64_t i = bw; i < (gnx - bw); i += 2)
            for (boost::uint64_t j = bw; j < (gnx - bw); j += 2)
                for (boost::uint64_t k = bw; k < (gnx - bw); k += 2)
                {
                    // Adjusted indices (for destination).
                    boost::uint64_t const id
                        = (bw + i) / 2 + idx.x() * ((gnx / 2) - bw);
                    boost::uint64_t const jd
                        = (bw + j) / 2 + idx.y() * ((gnx / 2) - bw);
                    boost::uint64_t const kd
                        = (bw + k) / 2 + idx.z() * ((gnx / 2) - bw);
   
                    if (derived().refine(e(id, jd, kd)))
                        return true;
                }
    
        return false;
    }

    bool unrefine(octree_server& e, child_index idx)
    {
        boost::uint64_t const bw = science().ghost_zone_width;
        boost::uint64_t const gnx = config().grid_node_length;

        for (boost::uint64_t i = bw; i < (gnx - bw); i += 2)
            for (boost::uint64_t j = bw; j < (gnx - bw); j += 2)
                for (boost::uint64_t k = bw; k < (gnx - bw); k += 2)
                {
                    // Adjusted indices (for destination).
                    boost::uint64_t const id
                        = (bw + i) / 2 + idx.x() * ((gnx / 2) - bw);
                    boost::uint64_t const jd
                        = (bw + j) / 2 + idx.y() * ((gnx / 2) - bw);
                    boost::uint64_t const kd
                        = (bw + k) / 2 + idx.z() * ((gnx / 2) - bw);

                    if (!derived().unrefine(e(id, jd, kd)))
                        return false;
                }

        return true;
    }

    refinement_criteria_base* clone() const
    {
        return new elementwise_refinement_criteria_base();
    }    

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & hpx::util::base_object_nonvirt<refinement_criteria_base>(*this);
    }
};

}

HPX_SERIALIZATION_REGISTER_TEMPLATE(
    (template <typename Derived>),
    (octopus::elementwise_refinement_criteria_base<Derived>)
)


#endif // OCTOPUS_EC4E7955_BF33_473F_9A44_FF274B50EEE8

