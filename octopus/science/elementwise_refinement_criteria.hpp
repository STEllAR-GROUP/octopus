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
struct elementwise_refinement_criteria_base
  : refinement_criteria_base
{
    elementwise_refinement_criteria_base()
    {
        hpx::actions::detail::guid_initialization<
            elementwise_refinement_criteria_base
        >();
    }

    static void register_base()
    {
        hpx::util::void_cast_register_nonvirt<
            elementwise_refinement_criteria_base, refinement_criteria_base
        >();
    }

    Derived& derived()
    {
        return static_cast<Derived&>(*this);
    }

    bool refine(octree_server& e, child_index idx)
    {
        boost::uint64_t const bw = science().ghost_zone_length;
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

                    array<double, 3> loc = e.center_coords(id, jd, kd);
   
                    if (derived().refine(e, e(id, jd, kd), loc))
                        return true;
                }
    
        return false;
    }

    bool unrefine(octree_server& e, child_index idx)
    {
        boost::uint64_t const bw = science().ghost_zone_length;
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

                    array<double, 3> loc = e.center_coords(id, jd, kd);

                    if (!derived().unrefine(e, e(id, jd, kd), loc))
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

