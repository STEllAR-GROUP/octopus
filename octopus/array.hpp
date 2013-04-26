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
#include <boost/move/move.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/serialization/access.hpp>

namespace octopus
{

template <typename T, typename Rep>
struct array_rep_traits 
{
    typedef T& reference;
    typedef T const& const_reference;

    template <typename Archive>
    static void serialize(
        Archive& ar
      , Rep& r
      , const unsigned int version
        )
    {
        ar & r; 
    }
};

template <typename T, boost::uint64_t Size>
struct array_rep_traits<T, boost::array<T, Size> > 
{
    typedef T& reference;
    typedef T const& const_reference;

    template <typename Archive>
    static void serialize(
        Archive& ar
      , boost::array<T, Size>& a
      , const unsigned int version
        )
    {
        ar & boost::serialization::make_array(a, Size);
    }
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
        array_rep_traits<T, Rep>::serialize(ar, rep_, version);
    }

  public:
    array() : data_()
    {
        for (size_type i = 0; i < Size; ++i)
            data_[i] = T();
    }

    array(Rep const& rep) : rep_(rep) {}

    array(BOOST_RV_REF(array) : rep_(boost::move(rhs.rep_)) {}

    array(array const& rhs) : rep_(rhs.rep) {}

    array(BOOST_RV_REF(array) rhs) : rep_(boost::move(rhs.rep_)) {}

    template <typename OtherRep>
    array& operator=(array<T, Size, OtherRep> const& rhs)
    {
        for (size_type i = 0; i < Size; ++i)
            rep_[i] = rhs_[i]; 
        return *this;
    }

    template <typename OtherRep>
    array& operator=(BOOST_RV_REF(array<T, Size, OtherRep>) rhs)
    {
        for (size_type i = 0; i < Size; ++i)
            rep_[i] = boost::move(rhs_[i]); 
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
    T const& s_;

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

    template <typename Archive>
    static void serialize(
        Archive& ar
      , array_dsel_scalar<T>& r
      , const unsigned int version
        )
    {
        ar & r; 
    }
};

#define OCTOPUS_DEFINE_ARRAY_ARITHMETIC_OPERATOR(OP, NAME)                    \
    template <typename T, typename LHS, typename RHS>                         \
    struct BOOST_PP_CAT(array_dsel_, NAME)                                    \
    {                                                                         \
        typedef boost::uint64_t size_type;                                    \
                                                                              \
      private:                                                                \
        LHS const& lhs_;                                                      \
        RHS const& rhs_;                                                      \
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
    template <typename T, boost::uint64_t Size, typename RHS>                 \
    array<T, Size, BOOST_PP_CAT(array_dsel_, NAME)                            \
        <T, array_dsel_scalar<T>, RHS> >                                      \
    operator OP (                                                             \
        T const& s                                                            \
      , array<T, Size, RHS> const& rhs                                        \
        )                                                                     \
    {                                                                         \
        typedef BOOST_PP_CAT(array_dsel_, NAME)                               \
            <T, array_dsel_scalar<T>, RHS> op;                                \
        return array<T, Size, op>(op(array_dsel_scalar<T>(s), rhs.rep()));    \
    }                                                                         \
                                                                              \
    template <typename T, boost::uint64_t Size, typename LHS>                 \
    array<T, Size, BOOST_PP_CAT(array_dsel_, NAME)                            \
        <T, array_dsel_scalar<T>, LHS> >                                      \
    operator OP (                                                             \
        array<T, Size, LHS> const& lhs                                        \
      , T const& s                                                            \
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
            lhs[i] = rhs[i];                                                  \
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
    Subject const& subject_;

  public:
    array_dsel_negation(Subject const& subject)
      : subject_(subject) 
    {}

    T operator[](size_type i) const
    {
        return -subject_[i]; 
    }
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

