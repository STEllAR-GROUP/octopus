////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/include/actions.hpp>
#include <hpx/lcos/future_wait.hpp>
#include <hpx/runtime/threads/thread_helpers.hpp>

#include <octopus/io/silo.hpp>
#include <octopus/math.hpp>
#include <octopus/engine/engine_interface.hpp>

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <cmath>

namespace octopus
{

void single_variable_silo_writer::start_write_locked(
    boost::uint64_t step
  , double time
  , std::string const& file
  , mutex_type::scoped_lock& l
    )
{
    OCTOPUS_ASSERT(l.owns_lock());

    // Make sure we closed the last epoch.
    stop_write_locked(l);

    step_ = step;
    time_ = time; 

    try
    {
        std::string s = boost::str( boost::format(file)
                                  % hpx::get_locality_id() % step_);
        file_ = DBCreate(s.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
    }
    // FIXME: Catch the specific boost.format exception.
    catch (...)
    {
        try
        {
            std::string s = boost::str( boost::format(file)
                                      % hpx::get_locality_id());
            file_ = DBCreate(s.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
        }
        // FIXME: Catch the specific boost.format exception.
        catch (...)
        {
            file_ = DBCreate(file.c_str()
                           , DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
        }
    }

    OCTOPUS_ASSERT(file_ != 0);

    directory_names_.reserve(config().levels_of_refinement + 1);
    for (boost::uint64_t i = 0; i < (config().levels_of_refinement + 1); ++i)
        directory_names_.push_back
            (boost::str(boost::format("/meshes_%1%") % i).c_str()); 

    for (boost::uint64_t i = 0; i < directory_names_.size(); ++i)
    {
        int error = DBMkDir(file_, directory_names_[i].c_str());
        OCTOPUS_ASSERT(error == 0);
    }
}

void single_variable_silo_writer::stop_write_locked(mutex_type::scoped_lock& l)
{
    OCTOPUS_ASSERT(l.owns_lock());

    merge_locked(l);

    if (file_)
    {
        DBClose(file_);
        file_ = 0; 
    }

    directory_names_.clear(); 
    merged_ = false;
}

void silo_perform_start_write(
    boost::uint64_t step
  , double time
  , std::string const& file
    )
{
    science().output.cast<single_variable_silo_writer>()->start_write
        (step, time, file);
}

void silo_perform_stop_write()
{
    science().output.cast<single_variable_silo_writer>()->stop_write();
}

}

HPX_PLAIN_ACTION(octopus::silo_perform_start_write
               , silo_perform_start_write_action);
HPX_ACTION_USES_MEDIUM_STACK(silo_perform_start_write_action);
HPX_PLAIN_ACTION(octopus::silo_perform_stop_write
               , silo_perform_stop_write_action);
HPX_ACTION_USES_MEDIUM_STACK(silo_perform_stop_write_action);

namespace octopus
{

void single_variable_silo_writer::begin_epoch(
    octree_server& e
  , double time
  , std::string const& file
    )
{
    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    silo_perform_start_write_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
    {
        if (file.empty())
            futures.emplace_back(hpx::async
                (act, targets[i], e.get_step(), time, file_name_));
        else
            futures.emplace_back(hpx::async
                (act, targets[i], e.get_step(), time, file));
    }

    hpx::wait(futures);
}

void single_variable_silo_writer::end_epoch(octree_server& e)
{
    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    silo_perform_stop_write_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
        futures.push_back(hpx::async(act, targets[i]));

    hpx::wait(futures);
}

void perform_merge(
    hpx::lcos::local::channel<void>& sync
  , DBfile* file 
  , std::vector<std::string> const& directory_names
  , std::string const& variable_name 
  , boost::uint64_t& step
  , double& time
    )
{ // {{{
    for (boost::uint64_t level = 0; level < directory_names.size(); ++level)
    {
        int error = DBSetDir(file, directory_names[level].c_str());
        OCTOPUS_ASSERT(error == 0);
    
        // Get a list of all the meshes and variables in that directory.
        DBtoc* contents = DBGetToc(file);

        // Make the mesh and variable names.
        boost::ptr_vector<char> mesh_names(contents->nqmesh);
        boost::ptr_vector<char> variable_names(contents->nqmesh);
        
        double grid_dim = config().spatial_domain;

        for (boost::uint64_t i = 1; i < level; ++i)
            grid_dim *= 0.5;

        for (boost::uint64_t j = 0; j < boost::uint64_t(contents->nqmesh); ++j)
        {
            std::string tmp;

            ///////////////////////////////////////////////////////////////////
            tmp  = directory_names[level];
            tmp += "/";
            tmp += contents->qmesh_names[j];

            // The extra character is for the terminating byte.
            mesh_names.push_back(new char[tmp.size() + 1]);
            std::strcpy(&mesh_names[j], tmp.c_str());

            ///////////////////////////////////////////////////////////////////
            tmp  = directory_names[level];
            tmp += "/";
            tmp += contents->qvar_names[j];

            // The extra character is for the terminating byte.
            variable_names.push_back(new char[tmp.size() + 1]);
            std::strcpy(&variable_names[j], tmp.c_str());
        }

        // when we change directories below, "contents" will go out of scope
        boost::uint64_t nqmesh = contents->nqmesh;

        error = DBSetDir(file, "/");
        OCTOPUS_ASSERT(error == 0);

        std::string multi_mesh_name
            = boost::str( boost::format("mesh_level_%1%")
                        % level);

        std::string multi_variable_name
            = boost::str( boost::format("%1%_level_%2%")
                        % variable_name % level); 


        {
            DBoptlist* optlist = DBMakeOptlist(4);
            DBObjectType type1 = DB_QUADRECT;
            DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type1);
            DBAddOption(optlist, DBOPT_CYCLE, &step);
            DBAddOption(optlist, DBOPT_DTIME, &time);

            error = DBPutMultimesh(file
                                 , multi_mesh_name.c_str()
                                 , nqmesh
                                 , mesh_names.c_array()
                                 , NULL, optlist);
            OCTOPUS_ASSERT(error == 0);
        }

        {
            DBoptlist* optlist = DBMakeOptlist(4);
            DBObjectType type1 = DB_QUADVAR;
            DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type1);
            DBAddOption(optlist, DBOPT_CYCLE, &step);
            DBAddOption(optlist, DBOPT_DTIME, &time);
            int type2 = DB_ROWMAJOR;
            DBAddOption(optlist, DBOPT_MAJORORDER, &type2);

            error = DBPutMultivar(file
                                , multi_variable_name.c_str()
                                , nqmesh
                                , variable_names.c_array()
                                , NULL, optlist);
            OCTOPUS_ASSERT(error == 0);
        }
    }

    sync.post();
} // }}}

void single_variable_silo_writer::merge_locked(mutex_type::scoped_lock& l)
{
    OCTOPUS_ASSERT(l.owns_lock());

    if (merged_ || !file_) return;

    hpx::lcos::local::channel<void> sync;
 
    hpx::applier::register_work_nullary(
        boost::bind(&perform_merge
                  , boost::ref(sync)
                  , file_
                  , boost::cref(directory_names_)
                  , boost::cref(variable_name_)
                  , boost::ref(step_)
                  , boost::ref(time_)
                    )
      , "perform_merge"
      , hpx::threads::pending
      , hpx::threads::thread_priority_normal
      , std::size_t(-1)
      , hpx::threads::thread_stacksize_medium
        );

    sync.get(); 
 
    merged_ = true;
}

void perform_write(
    hpx::lcos::local::channel<void>& sync
  , octree_server& e
  , DBfile* file 
  , std::vector<std::string> const& directory_names
  , std::string const& variable_name
  , boost::uint64_t variable_index
    )
{ // {{{
    boost::uint64_t const bw = science().ghost_zone_width;
    boost::uint64_t const gnx = config().grid_node_length;

    boost::uint64_t level = e.get_level();

    int nnodes[] = {
        int(gnx - 2 * bw + 1)
      , int(gnx - 2 * bw + 1)
      , int(gnx - 2 * bw + 1)
    };

    int nzones[] = {
        int(gnx - 2 * bw) 
      , int(gnx - 2 * bw) 
      , int(gnx - 2 * bw) 
    };

    //char* coordinate_names[] = { (char*) "X", (char*) "Y", (char*) "Z" };

    boost::scoped_array<double*> coordinates(new double* [3]);
    coordinates[0] = new double [nnodes[0]];
    coordinates[1] = new double [nnodes[1]];
    coordinates[2] = new double [nnodes[2]];

    boost::scoped_array<double> variables
        (new double [nzones[0] * nzones[1] * nzones[2]]);

    for (boost::uint64_t i = bw; i < (gnx - bw + 1); ++i) 
    {
        coordinates[0][i - bw] = e.x_face(i); 
        coordinates[1][i - bw] = e.y_face(i); 
        coordinates[2][i - bw] = e.z_face(i); 
    }

    for (boost::uint64_t i = bw; i < (gnx - bw); ++i) 
        for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
            for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
            {
                boost::uint64_t index = (i - bw)
                                      + (j - bw) * nzones[0]
                                      + (k - bw) * nzones[0] * nzones[1];
                variables[index] = e(i, j, k)[variable_index];
            }

    array<boost::uint64_t, 3> location = e.get_location();

    std::string mesh_name
        = boost::str( boost::format("mesh_L%i_%i_%i_%i")
                    % level
                    % location[0]
                    % location[1]
                    % location[2]);

    std::string value_name
        = boost::str( boost::format("%s_L%i_%i_%i_%i")
                    % variable_name
                    % level
                    % location[0]
                    % location[1]
                    % location[2]);

    int error = DBSetDir(file, directory_names[level].c_str());
    OCTOPUS_ASSERT(error == 0);

    {
        DBoptlist* optlist = DBMakeOptlist(1);
        // REVIEW: Verify this.
        int type = DB_ROWMAJOR;
        DBAddOption(optlist, DBOPT_MAJORORDER, &type);
        error = DBPutQuadmesh(file
                            , mesh_name.c_str()
                            //, coordinate_names
                            , NULL // SILO docs say this is ignored.
                            , coordinates.get()
                            , nnodes
                            , 3, DB_DOUBLE, DB_COLLINEAR, optlist);
        OCTOPUS_ASSERT(error == 0);
    }

    {
        DBoptlist* optlist = DBMakeOptlist(3);
        // REVIEW: Verify this.
        int type = DB_ROWMAJOR;
        DBAddOption(optlist, DBOPT_MAJORORDER, &type);
        error = DBPutQuadvar1(file
                            , value_name.c_str()
                            , mesh_name.c_str()
                            , variables.get()
                            , nzones
                            , 3, NULL, 0, DB_DOUBLE, DB_ZONECENT, optlist);
        OCTOPUS_ASSERT(error == 0);
    }

    delete[] coordinates[0];
    delete[] coordinates[1];
    delete[] coordinates[2]; 

    sync.post();
} // }}}

void single_variable_silo_writer::operator()(octree_server& e)
{
    mutex_type::scoped_lock l(mtx_);

    hpx::lcos::local::channel<void> sync;

    hpx::applier::register_work_nullary(
        boost::bind(&perform_write
                  , boost::ref(sync)
                  , boost::ref(e)
                  , file_
                  , boost::cref(directory_names_)
                  , boost::cref(variable_name_)
                  , boost::cref(variable_index_)
                    )
      , "perform_write"
      , hpx::threads::pending
      , hpx::threads::thread_priority_normal
      , std::size_t(-1)
      , hpx::threads::thread_stacksize_medium
        );

    sync.get(); 
}

}

