////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#define BOOST_NO_CXX11_SCOPED_ENUMS

#include <hpx/include/lcos.hpp>

#include <octopus/engine/engine_interface.hpp>
#include <octopus/filesystem.hpp>

namespace octopus
{

void call_here(hpx::util::function<void()> const& f) 
{
    f();  
}

}

HPX_PLAIN_ACTION(octopus::call_here, call_here_action);

namespace octopus
{

std::vector<hpx::future<void> > call_everywhere(
    hpx::util::function<void()> const& f
    ) 
{
    OCTOPUS_ASSERT_MSG(!localities().empty(),
                       "no localities supporting Octopus available");

    std::vector<hpx::future<void> > calls;
    calls.reserve(localities().size());

    for (boost::uint64_t i = 0; i < localities().size(); ++i)
        calls.emplace_back(hpx::async<call_here_action>(localities()[i], f));

    return calls;
}

struct backup_checkpoint_locally
{
  private:
    std::string suffix_;

  public:
    backup_checkpoint_locally() : suffix_() {}

    backup_checkpoint_locally(std::string const& suffix)
      : suffix_(suffix)
    {}

    void operator()() const
    {
        checkpoint().flush();

        std::string file_name;

        try
        {
            file_name = boost::str( boost::format(config().checkpoint_file)
                                  % hpx::get_locality_id());
        }
        // FIXME: Catch the specific boost.format exception.
        catch (...)
        {
            file_name = config().checkpoint_file;
        }

        namespace fs = boost::filesystem;

        fs::path from(file_name);

        fs::path to(file_name + suffix_);

        fs::copy_file(from, to, fs::copy_option::overwrite_if_exists);
    }

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int)
    {
        ar & suffix_;
    };
};

void backup_checkpoint(
    std::string const& suffix
    )
{
    backup_checkpoint_locally bcl(suffix);
    hpx::wait(call_everywhere(bcl));
}


}


