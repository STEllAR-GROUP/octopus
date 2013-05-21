////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_2B645EED_85D3_4967_89C9_0905C684BE8C)
#define OCTOPUS_2B645EED_85D3_4967_89C9_0905C684BE8C

#include <octopus/assert.hpp>

#include <boost/array.hpp>
#include <boost/cstdint.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/move/move.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>

// TODO: Split array_rep_traits into two traits classes..

namespace octopus
{

template <typename T, typename Rep>
struct array_rep_traits 
{
    typedef T& reference;
    typedef T const& const_reference;

    typedef Rep const& storage_type;
};

template <
    typename T
  , boost::uint64_t Size
  , typename Rep = boost::array<T, Size>
>
struct array
{
    typedef boost::uint64_t size_type;

  private:
    Rep rep_;

    BOOST_COPYABLE_AND_MOVABLE(array);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & rep_;
    }

  public:
    array() : rep_()
    {
        for (size_type i = 0; i < Size; ++i)
            rep_[i] = T();
    }

    array(Rep const& rep) : rep_(rep) {}

    array(BOOST_RV_REF(Rep) rep) : rep_(boost::move(rep)) {}

    array(array const& rhs) : rep_(rhs.rep_) {}

    array(BOOST_RV_REF(array) rhs) : rep_(boost::move(rhs.rep_)) {}
 
    template <typename OtherRep>
    array(array<T, Size, OtherRep> const& rhs)
    {
        *this = rhs;
    }

    template <typename OtherRep>
    array(BOOST_RV_REF_3_TEMPL_ARGS(array, T, Size, OtherRep) rhs)
    {
        *this = boost::move(rhs);
    }

    array& operator=(BOOST_COPY_ASSIGN_REF(array) rhs)
    {
        rep_ = rhs.rep_;
        return *this;
    }

    array& operator=(BOOST_RV_REF(array) rhs)
    {
        rep_ = boost::move(rhs.rep_);
        return *this;
    }

    template <typename OtherRep>
    array& operator=(
        BOOST_COPY_ASSIGN_REF_3_TEMPL_ARGS(array, T, Size, OtherRep) rhs
        )
    {
        for (size_type i = 0; i < Size; ++i)
            rep_[i] = rhs[i]; 
        return *this;
    }

    template <typename OtherRep>
    array& operator=(
        BOOST_RV_REF_3_TEMPL_ARGS(array, T, Size, OtherRep) rhs
        )
    {
        for (size_type i = 0; i < Size; ++i)
            rep_[i] = boost::move(rhs[i]); 
        return *this;
    }

    array& operator=(T const& value)
    {
        for (size_type i = 0; i < Size; ++i)
            rep_[i] = value;
        return *this;
    }

    size_type size() const
    {
        return Size;
    } 

    Rep const& rep() const
    {
        return rep_;
    } 
    
    typename array_rep_traits<T, Rep>::reference
    operator()(size_type c)
    {
        OCTOPUS_ASSERT_FMT_MSG(c < Size,
            "coordinate (%1%) is larger than the length (%2%)",
            c % Size);  
        return rep_[c];
    }

    typename array_rep_traits<T, Rep>::const_reference
    operator()(size_type c) const
    {
        OCTOPUS_ASSERT_FMT_MSG(c < Size,
            "coordinate (%1%) is larger than the length (%2%)",
            c % Size);  
        return rep_[c];
    }

    typename array_rep_traits<T, Rep>::reference
    operator[](size_type c)
    {
        return (*this)(c);
    }

    typename array_rep_traits<T, Rep>::const_reference
    operator[](size_type c) const
    {
        return (*this)(c);
    }
};

template <typename T>
struct array_dsel_scalar
{
    typedef boost::uint64_t size_type;

  private:
    T s_;

  public:
    array_dsel_scalar(T const& s) : s_(s) {}

    T const& operator[] (size_type) const
    {
        return s_;
    }
};

template <typename T>
struct array_rep_traits<T, array_dsel_scalar<T> >  
{
    typedef T const& reference;
    typedef T const& const_reference;

    typedef array_dsel_scalar<T> storage_type;
};

#define OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR(OP, NAME)                    \
    template <typename T, typename LHS, typename RHS>                         \
    struct BOOST_PP_CAT(array_dsel_, NAME)                                    \
    {                                                                         \
        typedef boost::uint64_t size_type;                                    \
                                                                              \
      private:                                                                \
        typename array_rep_traits<T, LHS>::storage_type lhs_;                 \
        typename array_rep_traits<T, RHS>::storage_type rhs_;                 \
                                                                              \
      public:                                                                 \
        BOOST_PP_CAT(array_dsel_, NAME)(LHS const& lhs, RHS const& rhs)       \
          : lhs_(lhs), rhs_(rhs)                                              \
        {}                                                                    \
                                                                              \
        T operator[](size_type i) const                                       \
        {                                                                     \
            return lhs_[i] OP rhs_[i];                                        \
        }                                                                     \
    };                                                                        \
                                                                              \
    template <typename T, typename LHS, typename RHS>                         \
    struct array_rep_traits<T, BOOST_PP_CAT(array_dsel_, NAME)<T, LHS, RHS> > \
    {                                                                         \
        typedef T reference;                                                  \
        typedef T const_reference;                                            \
                                                                              \
        typedef BOOST_PP_CAT(array_dsel_, NAME)<T, LHS, RHS> const&           \
            storage_type;                                                     \
    };                                                                        \
                                                                              \
    template <typename T, boost::uint64_t Size, typename LHS, typename RHS>   \
    array<T, Size, BOOST_PP_CAT(array_dsel_, NAME)<T, LHS, RHS> >             \
    operator OP (                                                             \
        array<T, Size, LHS> const& lhs                                        \
      , array<T, Size, RHS> const& rhs                                        \
        )                                                                     \
    {                                                                         \
        typedef BOOST_PP_CAT(array_dsel_, NAME)<T, LHS, RHS> op;              \
        return array<T, Size, op>(op(lhs.rep(), rhs.rep()));                  \
    }                                                                         \
                                                                              \
    template <typename I, typename T, boost::uint64_t Size, typename RHS>     \
    typename boost::enable_if_c<boost::is_arithmetic<I>::value,               \
        array<T, Size, BOOST_PP_CAT(array_dsel_, NAME)                        \
            <T, array_dsel_scalar<T>, RHS> >                                  \
    >::type                                                                   \
    operator OP (                                                             \
        I s                                                                   \
      , array<T, Size, RHS> const& rhs                                        \
        )                                                                     \
    {                                                                         \
        typedef BOOST_PP_CAT(array_dsel_, NAME)                               \
            <T, array_dsel_scalar<T>, RHS> op;                                \
        return array<T, Size, op>(op(array_dsel_scalar<T>(s), rhs.rep()));    \
    }                                                                         \
                                                                              \
    template <typename I, typename T, boost::uint64_t Size, typename LHS>     \
    typename boost::enable_if_c<boost::is_arithmetic<I>::value,               \
        array<T, Size, BOOST_PP_CAT(array_dsel_, NAME)                        \
            <T, LHS, array_dsel_scalar<T> > >                                 \
    >::type                                                                   \
    operator OP (                                                             \
        array<T, Size, LHS> const& lhs                                        \
      , I s                                                                   \
        )                                                                     \
    {                                                                         \
        typedef BOOST_PP_CAT(array_dsel_, NAME)                               \
            <T, LHS, array_dsel_scalar<T> > op;                               \
        return array<T, Size, op>(op(lhs.rep(), array_dsel_scalar<T>(s)));    \
    }                                                                         \
                                                                              \
    template <typename T, boost::uint64_t Size, typename Rep>                 \
    array<T, Size>                                                            \
    operator BOOST_PP_CAT(OP, =) (                                            \
        array<T, Size>& lhs                                                   \
      , array<T, Size, Rep> const& rhs                                        \
        )                                                                     \
    {                                                                         \
        for (typename array<T, Size>::size_type i = 0; i < lhs.size(); ++i)   \
            lhs[i] BOOST_PP_CAT(OP, =) rhs[i];                                \
        return lhs;                                                           \
    }                                                                         \
                                                                              \
    template <typename I, typename T, boost::uint64_t Size>                   \
    typename boost::enable_if_c<boost::is_arithmetic<I>::value,               \
        array<T, Size>                                                        \
    >::type                                                                   \
    operator BOOST_PP_CAT(OP, =) (                                            \
        array<T, Size>& lhs                                                   \
      , I rhs                                                                 \
        )                                                                     \
    {                                                                         \
        for (typename array<T, Size>::size_type i = 0; i < lhs.size(); ++i)   \
            lhs[i] BOOST_PP_CAT(OP, =) rhs;                                   \
        return lhs;                                                           \
    }                                                                         \
    /**/

OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR(-, minus)
OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR(+, plus)
OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR(*, multiply)
OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR(/, divides)

#undef OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR

template <typename T, typename Subject>
struct array_dsel_negation 
{
    typedef boost::uint64_t size_type;

  private:
    typename array_rep_traits<T, Subject>::storage_type subject_;

  public:
    array_dsel_negation(Subject const& subject)
      : subject_(subject) 
    {}

    T operator[](size_type i) const
    {
        return -subject_[i]; 
    }
};

template <typename T, typename Subject>
struct array_rep_traits<T, array_dsel_negation<T, Subject> >  
{
    typedef T reference;
    typedef T const_reference;

    typedef array_dsel_negation<T, Subject> const& storage_type;
};

template <typename T, boost::uint64_t Size, typename Subject>
array<T, Size, Subject> operator-(
    array<T, Size, Subject> const& subject
    ) 
{
    typedef array_dsel_negation<T, Subject> op;
    return array<T, Size, op>(op(subject.rep()));
}

}

#endif // OCTOPUS_2B645EED_85D3_4967_89C9_0905C684BE8C

