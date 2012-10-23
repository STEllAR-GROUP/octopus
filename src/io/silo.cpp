////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Bryce Adelstein-Lelbach
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <octopus/io/silo.hpp>
#include <octopus/engine/engine_interface.hpp>

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

namespace octopus
{

void single_variable_silo_writer::open_locked(
    boost::uint64_t step
  , mutex_type::scoped_lock& l
    )
{
    OCTOPUS_ASSERT(l.owns_lock());

    merge_locked(l);
    close_locked(l);

    step_ = step;

    file_ = DBCreate(boost::str( boost::format(file_name_)
                               % step_ % hpx::get_locality_id()).c_str() 
                   , DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
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

void single_variable_silo_writer::close_locked(mutex_type::scoped_lock& l)
{
    OCTOPUS_ASSERT(l.owns_lock());

    if (file_)
    {
        DBClose(file_);
        file_ = 0; 
    }

    directory_names_.clear(); 
}

void single_variable_silo_writer::merge_locked(mutex_type::scoped_lock& l)
{
    OCTOPUS_ASSERT(l.owns_lock());

    if (merged_ || !file_) return;
  
    for (boost::uint64_t level = 0; level <= directory_names_.size(); ++level)
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

        DBSetDir(file_, "/");

        std::string multi_mesh_name
            = boost::str( boost::format("mesh_level_%1%")
                        % level);

        std::string multi_variable_name
            = boost::str( boost::format("%1%_level_%2%")
                        % variable_name_ % level); 


        {
            DBoptlist* optlist = DBMakeOptlist(1);
            DBObjectType type = DB_QUADRECT;
            DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type);

            int error = DBPutMultimesh(file_
                                     , multi_mesh_name.c_str()
                                     , contents->nqmesh
                                     , mesh_names.c_array()
                                     , NULL, optlist);
            OCTOPUS_ASSERT(error == 0);
        }

        {
            DBoptlist* optlist = DBMakeOptlist(1);
            DBObjectType type = DB_QUADVAR;
            DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type);

            int error = DBPutMultivar(file_
                                    , multi_variable_name.c_str()
                                    , contents->nqmesh
                                    , variable_names.c_array()
                                    , NULL, optlist);
            OCTOPUS_ASSERT(error == 0);
        }
    }

    merged_ = true;
}

void single_variable_silo_writer::operator()(octree_server& e)
{
    mutex_type::scoped_lock(mtx_);

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

    char* coordinate_names[] = { (char*) "X", (char*) "Y", (char*) "Z" };

    boost::ptr_vector<double> coordinates(nnodes[0] * nnodes[0] * nnodes[0]);
    boost::scoped_array<double*> variables
        (new double* [nzones[0] * nzones[0] * nzones[0]]);

    for (boost::uint64_t i = bw; i < (gnx - bw + 1); ++i) 
        for (boost::uint64_t j = bw; j < (gnx - bw + 1); ++j) 
            for (boost::uint64_t k = bw; k < (gnx - bw + 1); ++k) 
            {
                coordinates.push_back(new double[3]);

                (&coordinates.back())[0] = e.xf(i); 
                (&coordinates.back())[1] = e.yf(i); 
                (&coordinates.back())[2] = e.zf(i); 
            }

    for (boost::uint64_t i = bw; i < (gnx - bw); ++i) 
        for (boost::uint64_t j = bw; j < (gnx - bw); ++j) 
            for (boost::uint64_t k = bw; k < (gnx - bw); ++k) 
            {
                boost::uint64_t index = (i - bw)
                                      + (j - bw) * nzones[0]
                                      + (k - bw) * nzones[0] * nzones[1];
                variables[index] = &e(i, j, k)[variable_index_];
            }

    std::string mesh_name
        = boost::str( boost::format("mesh_%1%_%2%")
                    % e.get_gid().get_msb()
                    % e.get_gid().get_lsb());

    std::string variable_name
        = boost::str( boost::format("%1%_%2%_%3%")
                    % variable_name_
                    % e.get_gid().get_msb()
                    % e.get_gid().get_lsb());

    int error = DBSetDir(file_, directory_names_[level].c_str());
    OCTOPUS_ASSERT(error == 0);

    error = DBPutQuadmesh(file_
                        , mesh_name.c_str()
                        , coordinate_names
                        , coordinates.c_array()
                        , nnodes
                        , 3, DB_DOUBLE, DB_COLLINEAR, NULL);
    OCTOPUS_ASSERT(error == 0);

    error = DBPutQuadvar1(file_
                        , variable_name.c_str()
                        , mesh_name.c_str()
                        , variables.get()
                        , nzones
                        , 3, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
    OCTOPUS_ASSERT(error == 0);
}

}

