////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2007-2012 Hartmut Kaiser
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#if !defined(OCTOPUS_630F2401_C82B_4EBA_B213_2FD6E303E050)
#define OCTOPUS_630F2401_C82B_4EBA_B213_2FD6E303E050

#include <hpx/lcos/detail/future_data.hpp>
#include <hpx/lcos/local/packaged_continuation.hpp>

namespace octopus
{

/// An asynchronous, single value channel
template <typename T>
struct channel
{
  private:
    typedef hpx::lcos::detail::future_data<T> future_data;

    boost::intrusive_ptr<future_data> data_;

    BOOST_COPYABLE_AND_MOVABLE(vector3d);

  public:
    typedef typename future_data::completed_callback_type
        completed_callback_type;

    channel() : data_(new future_data())
    {}

    channel(channel const& other) : data_()
    {
        *this = other;
    }

    channel(BOOST_RV_REF(channel) other) : data_()
    {
        *this = boost::move(other);
    }

    explicit channel(BOOST_RV_REF(T) init) : data_(new future_data())
    {
        data_->set_data(init); 
    }

    explicit channel(T const& init) : data_(new future_data())
    {
        data_->set_data(init); 
    }

    ~channel()
    {
        if (data_)
        {
            if (data_->is_ready())
            {
                data_->set_error(hpx::broken_promise,
                    "channel<T>::~channel",
                    "deleting owner before channel value has been consumed");
            }

            data_->deleting_owner();
        }
    }

    channel& operator=(BOOST_COPY_ASSIGN_REF(channel) other)
    {
        OCTOPUS_ASSERT(data_);

        if (this != &other)
        {
/*
            if (data_->is_ready())
            {
                data_->set_error(hpx::broken_promise,
                    "channel<T>::operator=()",
                    "deleting owner before channel value has been consumed");
            }

            data_->deleting_owner();
*/
            data_ = other.data_;
        }

        return *this;
    }

    channel& operator=(BOOST_RV_REF(channel) other)
    {
        OCTOPUS_ASSERT(data_);

        if (this != &other)
        {
/*
            if (data_->is_ready())
            {
                data_->set_error(hpx::broken_promise,
                    "channel<T>::operator=()",
                    "deleting owner before channel value has been consumed");
            }

            data_->deleting_owner();
*/
            data_ = boost::move(other.data_);
            other.data_.reset();
        }

        return *this;
    }

    void swap(channel& other)
    {
        data_.swap(other.data_);
    }

    void reset()
    {
        OCTOPUS_ASSERT(data_);

/*
        if (data_->is_ready())
        {
            data_->set_error(hpx::broken_promise,
                "channel<T>::clear()",
                "clearing owner before channel value has been retrieved");
        }

        data_->deleting_owner();
*/ 

        data_->reset();
   }

    T take(hpx::error_code& ec = hpx::throws) 
    {
        OCTOPUS_ASSERT(data_);
        OCTOPUS_ASSERT(data_->ready());
        T tmp = data_->move_data(ec);
        data_->reset();
        return boost::move(tmp);
    }

    void post(BOOST_RV_REF(T) result)
    {
        OCTOPUS_ASSERT(data_);
        data_->set_data(result);
    }

    void post(T const& result)
    {
        OCTOPUS_ASSERT(data_);
        data_->set_data(result);
    }

    template <typename F>
    hpx::future<typename boost::result_of<F(hpx::future<T>)>::type>
    then(BOOST_FWD_REF(F) f)
    {
        /*
        typedef typename boost::result_of<F(hpx::future<T>)>::type result_type;

        // Create a continuation.
        typedef hpx::lcos::local::packaged_continuation<result_type, T, T>
            cont_type;

        boost::shared_ptr<cont_type> p(
            boost::make_shared<cont_type>(
                hpx::util::bind(boost::forward<F>(f)
                              , hpx::util::placeholders::_1)
            )
        );

        // Bind a on_completed handler to the future data which will invoke the
        // continuation.
        data_->set_on_completed(
            hpx::util::bind(&cont_type::on_value_ready
                          , p, hpx::util::placeholders::_1));

        // The continuation will keep the future alive. We'll need a new one.
        //data_ = new future_data();

        return p->get_future();
        */
        return hpx::future<T>(data_).when
            (boost::forward<completed_callback_type>(f));
    }

    bool ready() const
    {
        OCTOPUS_ASSERT(data_);
        return data_->ready();
    }
};

}

#endif // OCTOPUS_630F2401_C82B_4EBA_B213_2FD6E303E050

