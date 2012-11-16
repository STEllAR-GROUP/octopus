////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Zach Byerly
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/include/actions.hpp>
#include <hpx/lcos/future_wait.hpp>

#include <octopus/io/silo.hpp>
#include <octopus/engine/engine_interface.hpp>

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

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
        file_ = DBCreate(boost::str( boost::format(file)
                                   % hpx::get_locality_id() % step_).c_str() 
                       , DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
    }
    // FIXME: Catch the specific boost.format exception.
    catch (...)
    {
        try
        {
            file_ = DBCreate(boost::str( boost::format(file)
                                       % hpx::get_locality_id()).c_str() 
                           , DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
        }
        // FIXME: Catch the specific boost.format exception.
        catch (...)
        {
            file_ = DBCreate(file.c_str()
                           , DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
        }
    }
 
    OCTOPUS_ASSERT(file_ != 0);

    directory_names_.reserve(config().max_refinement_level + 1);
    for (boost::uint64_t i = 0; i < (config().max_refinement_level + 1); ++i)
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

void start_write_locally(
    boost::uint64_t step
  , double time
  , std::string const& file
    )
{
    science().output.cast<single_variable_silo_writer>()->start_write
        (step, time, file);
}

void stop_write_locally()
{
    science().output.cast<single_variable_silo_writer>()->stop_write();
}

}

HPX_PLAIN_ACTION(octopus::start_write_locally, start_write_locally_action);
HPX_PLAIN_ACTION(octopus::stop_write_locally, stop_write_locally_action);

namespace octopus
{

void single_variable_silo_writer::begin_epoch(
    octree_server& e
  , std::string const& file
    )
{
    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    start_write_locally_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
        futures.push_back(hpx::async
            (act, targets[i], e.get_step(), e.get_time(), file));

    hpx::wait(futures);
}

void single_variable_silo_writer::end_epoch(octree_server& e)
{
    std::vector<hpx::id_type> const& targets = localities();

    std::vector<hpx::future<void> > futures;
    futures.reserve(targets.size());

    stop_write_locally_action act;

    for (boost::uint64_t i = 0; i < targets.size(); ++i)
        futures.push_back(hpx::async(act, targets[i]));

    hpx::wait(futures);
}

void single_variable_silo_writer::merge_locked(mutex_type::scoped_lock& l)
{
    OCTOPUS_ASSERT(l.owns_lock());

    if (merged_ || !file_) return;
  
    for (boost::uint64_t level = 0; level < directory_names_.size(); ++level)
    {
        int error = DBSetDir(file_, directory_names_[level].c_str());
        OCTOPUS_ASSERT(error == 0);
    
        // Get a list of all the meshes and variables in that directory.
        DBtoc* contents = DBGetToc(file_);

        // Make the mesh and variable names.
        boost::ptr_vector<char> mesh_names(contents->nqmesh);
        boost::ptr_vector<char> variable_names(contents->nqmesh);
        
        for (boost::uint64_t j = 0; j < boost::uint64_t(contents->nqmesh); ++j)
        {
            std::string tmp;

            ///////////////////////////////////////////////////////////////////
            tmp  = directory_names_[level];
            tmp += "/";
            tmp += contents->qmesh_names[j];

            // The extra character is for the terminating byte.
            mesh_names.push_back(new char[tmp.size() + 1]);
            std::strcpy(&mesh_names[j], tmp.c_str());

            ///////////////////////////////////////////////////////////////////
            tmp  = directory_names_[level];
            tmp += "/";
            tmp += contents->qvar_names[j];

            // The extra character is for the terminating byte.
            variable_names.push_back(new char[tmp.size() + 1]);
            std::strcpy(&variable_names[j], tmp.c_str());
        }

        // when we change directories below, "contents" will go out of scope
        boost::uint64_t nqmesh = contents->nqmesh;

        error = DBSetDir(file_, "/");
        OCTOPUS_ASSERT(error == 0);

        std::string multi_mesh_name
            = boost::str( boost::format("mesh_level_%1%")
                        % level);

        std::string multi_variable_name
            = boost::str( boost::format("%1%_level_%2%")
                        % variable_name_ % level); 


        {
            DBoptlist* optlist = DBMakeOptlist(4);
            DBObjectType type1 = DB_QUADRECT;
            DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type1);
            DBAddOption(optlist, DBOPT_CYCLE, &step_);
            DBAddOption(optlist, DBOPT_DTIME, &time_);

            error = DBPutMultimesh(file_
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
            DBAddOption(optlist, DBOPT_CYCLE, &step_);
            DBAddOption(optlist, DBOPT_DTIME, &time_);
            int type2 = DB_ROWMAJOR;
            DBAddOption(optlist, DBOPT_MAJORORDER, &type2);

            error = DBPutMultivar(file_
                                , multi_variable_name.c_str()
                                , nqmesh
                                , variable_names.c_array()
                                , NULL, optlist);
            OCTOPUS_ASSERT(error == 0);
        }
    }

    merged_ = true;
}

void single_variable_silo_writer::operator()(octree_server& e)
{
    mutex_type::scoped_lock l(mtx_);

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
                variables[index] = e(i, j, k)[variable_index_];
            }

    std::string mesh_name
        = boost::str( boost::format("mesh_%016x%016x")
                    % e.get_gid().get_msb()
                    % e.get_gid().get_lsb());

    std::string variable_name
        = boost::str( boost::format("%s_%016x%016x")
                    % variable_name_
                    % e.get_gid().get_msb()
                    % e.get_gid().get_lsb());

    int error = DBSetDir(file_, directory_names_[level].c_str());
    OCTOPUS_ASSERT(error == 0);

    {
        DBoptlist* optlist = DBMakeOptlist(1);
        // REVIEW: Verify this.
        int type = DB_ROWMAJOR;
        DBAddOption(optlist, DBOPT_MAJORORDER, &type);
        error = DBPutQuadmesh(file_
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
        error = DBPutQuadvar1(file_
                            , variable_name.c_str()
                            , mesh_name.c_str()
                            , variables.get()
                            , nzones
                            , 3, NULL, 0, DB_DOUBLE, DB_ZONECENT, optlist);
        OCTOPUS_ASSERT(error == 0);
    }

    delete[] coordinates[0];
    delete[] coordinates[1];
    delete[] coordinates[2]; 
}

}

