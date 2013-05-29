////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_EEC57C22_4221_4E32_8288_DC2A4F9D957B)
#define OCTOPUS_EEC57C22_4221_4E32_8288_DC2A4F9D957B

#include <octopus/traits.hpp>

#include <hpx/util/static.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/lcos/future_wait.hpp>

namespace octopus
{

template <typename T, typename Tag>
void update_here(typename parameter_type<T>::type t)
{
    hpx::util::static_<T, Tag> storage;
    storage.get() = t;
}
    
template <typename T, typename Tag>
struct update_here_action
  : hpx::actions::make_action<
        void (*)(typename parameter_type<T>::type)
      , &update_here<T, Tag>
      , update_here_action<T, Tag>
    >
{};

template <typename T, typename Tag>
struct global_variable
{
  private:
    void update_everywhere(typename parameter_type<T>::type t)
    {
        std::vector<hpx::id_type> targets = hpx::find_all_localities();

        std::vector<hpx::future<void> > futures;
        futures.reserve(targets.size()); 

        for (boost::uint64_t i = 0; i < targets.size(); ++i)
        {
            futures.push_back(hpx::async<update_here_action<T, Tag> >
                (targets[i], t)); 
        }

        hpx::lcos::wait(futures); 
    }

  public:
    global_variable() {}

    global_variable(typename parameter_type<T>::type t)
    {
        hpx::util::static_<T, Tag> storage;
        storage.get() = t; 
    }

    global_variable& operator=(typename parameter_type<T>::type t)
    {
        update_everywhere(t);
        return *this;
    }

    operator T const&() const
    {
        hpx::util::static_<T, Tag> storage; 
        return storage.get();
    }

    T const& get() const
    {
        hpx::util::static_<T, Tag> storage; 
        return storage.get();
    }
};

/// Let the INI config reader get to the real type.
template <typename T, typename Tag>
struct proxied_type<global_variable<T, Tag> >
{
    typedef T type;
};

/// Ensure that the INI config reader handles bools properly.
template <typename Tag>
struct is_bool<global_variable<bool, Tag> > : boost::mpl::true_ {};


/// The compiler won't do the right thing in this (e.g. call the conversion
/// operator). Necessary for the INI config reader.
template <typename T, typename Tag>
inline std::ostream& operator<<(
    std::ostream& os
  , global_variable<T, Tag> gv
    )
{
    return os << gv.get();
}

}

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template <typename T, typename Tag>),
    (octopus::update_here_action<T, Tag>))

/// Usage: OCTOPUS_GLOBAL_VARIABLE((type), name)
/// Usage: OCTOPUS_GLOBAL_VARIABLE((type), name, (initial_value))
#define OCTOPUS_GLOBAL_VARIABLE(...)                                          \
    HPX_UTIL_EXPAND_(BOOST_PP_CAT(                                            \
        OCTOPUS_GLOBAL_VARIABLE_, HPX_UTIL_PP_NARG(__VA_ARGS__)               \
    )(__VA_ARGS__))                                                           \
    /**/

#define OCTOPUS_GLOBAL_VARIABLE_2(type, name)                                 \
struct BOOST_PP_CAT(                                                          \
    BOOST_PP_CAT(__hpx_global_variable_tag_, __LINE__),                       \
    BOOST_PP_CAT(_, name)) {};                                                \
                                                                              \
octopus::global_variable<HPX_UTIL_STRIP(type), BOOST_PP_CAT(                  \
    BOOST_PP_CAT(__hpx_global_variable_tag_, __LINE__),                       \
    BOOST_PP_CAT(_, name))> name =                                            \
        octopus::global_variable<HPX_UTIL_STRIP(type), BOOST_PP_CAT(          \
            BOOST_PP_CAT(__hpx_global_variable_tag_, __LINE__),               \
            BOOST_PP_CAT(_, name))>();                                        \
    /**/

#define OCTOPUS_GLOBAL_VARIABLE_3(type, name, initial_value)                  \
struct BOOST_PP_CAT(                                                          \
    BOOST_PP_CAT(__hpx_global_variable_tag_, __LINE__),                       \
    BOOST_PP_CAT(_, name)) {};                                                \
                                                                              \
octopus::global_variable<HPX_UTIL_STRIP(type), BOOST_PP_CAT(                  \
    BOOST_PP_CAT(__hpx_global_variable_tag_, __LINE__),                       \
    BOOST_PP_CAT(_, name))> name =                                            \
        octopus::global_variable<HPX_UTIL_STRIP(type), BOOST_PP_CAT(          \
            BOOST_PP_CAT(__hpx_global_variable_tag_, __LINE__),               \
            BOOST_PP_CAT(_, name))>(HPX_UTIL_STRIP(initial_value));           \
    /**/

#endif // OCTOPUS_EEC57C22_4221_4E32_8288_DC2A4F9D957B

