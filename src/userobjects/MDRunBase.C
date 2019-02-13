/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MDRunBase.h"
#include "MooseMesh.h"

#include "libmesh/mesh_tools.h"

// custom data load and data store methods for MDParticles class
template <>
inline void
dataStore(std::ostream & stream, MDRunBase::MDParticles & pl, void * context)
{
  dataStore(stream, pl.pos, context);
  dataStore(stream, pl.vel, context);
  dataStore(stream, pl.id, context);
  dataStore(stream, pl.elem_id, context);
}

template <>
inline void
dataLoad(std::istream & stream, MDRunBase::MDParticles & pl, void * context)
{
  dataLoad(stream, pl.pos, context);
  dataLoad(stream, pl.vel, context);
  dataLoad(stream, pl.id, context);
  dataLoad(stream, pl.elem_id, context);
}

template <>
InputParameters
validParams<MDRunBase>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  params.suppressParameter<ExecFlagEnum>("execute_on");

  params.addClassDescription(
      "Base class for execution of coupled molecular dynamics MOOSE calculations.");
  return params;
}

MDRunBase::MDRunBase(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension()),
    _nproc(_app.n_processors()),
    _bbox(_nproc)
{
}

void
MDRunBase::initialSetup()
{
  for (unsigned int j = 0; j < _nproc; ++j)
    _bbox[j] = MeshTools::create_processor_bounding_box(_mesh, j);
}

void
MDRunBase::timestepSetup()
{
  // update/init subdomain bounding boxes
  for (unsigned int j = 0; j < _nproc; ++j)
    _bbox[j] = MeshTools::create_processor_bounding_box(_mesh, j);

  // callback for updating md particle list
  updateParticleList();
}

void
MDRunBase::updateKDTree()
{
  _kd_tree = libmesh_make_unique<KDTree>(_md_particles.pos, 50);
}

const std::vector<unsigned int>
MDRunBase::elemParticles(unique_id_type elem_id) const
{
  if (_elem_particles.find(elem_id) != _elem_particles.end())
    return _elem_particles.find(elem_id)->second;
  return std::vector<unsigned int>(0);
}

void
MDRunBase::mapMDParticles()
{
  // clear data structures
  _elem_particles.clear();
  _md_particles.elem_id.assign(_md_particles.pos.size(), libMesh::DofObject::invalid_unique_id);

  // loop over semi-local elements to ensure consistent handling of
  // points on processor boundaries
  const libMesh::MeshBase & mesh_base = _mesh.getMesh();
  for (const auto & elem : as_range(mesh_base.active_semilocal_elements_begin(),
                                    mesh_base.active_semilocal_elements_end()))
  {
    // find all points within an inflated bounding box
    std::vector<std::pair<std::size_t, Real>> indices_dist;
    BoundingBox bbox = elem->loose_bounding_box();
    Point center = 0.5 * (bbox.min() + bbox.max());
    Real radius = (bbox.max() - center).norm();
    _kd_tree->radiusSearch(center, radius, indices_dist);

    for (auto & p : indices_dist)
    {
      Point candidate = _md_particles.pos[p.first];

      // avoid double counting of elements on element boundaries, smallest
      // elem id gets to own the point, this is consistent across processors
      // since we loop over semi-local elements
      if (_md_particles.elem_id[p.first] != libMesh::DofObject::invalid_unique_id &&
          elem->unique_id() > _md_particles.elem_id[p.first])
        continue;

      // contains_point performs cheap bounding box test, hence no need to do it before,
      // serious candidates need to go through the expensive test
      if (elem->contains_point(candidate))
      {
        // convenience variable for the MD particle we deal with
        unsigned int pp = p.first;

        // remove this particle from the entry for the old element to avoid double
        // counting
        if (_md_particles.elem_id[pp] != libMesh::DofObject::invalid_unique_id)
        {
          // get the stored unique id of the _md_particle
          unique_id_type old_elem_id = _md_particles.elem_id[pp];

          // entry in _elem_particles should exist but guard w/ assert
          mooseAssert(_elem_particles.find(old_elem_id) != _elem_particles.end(),
                      "Entry " << old_elem_id << " in _elem_particles should exist.");
          std::vector<unsigned int> id_list = _elem_particles.find(old_elem_id)->second;

          // go through this list and find the entry with the right MD particle entry
          // and remove it
          for (unsigned int i = 0; i < id_list.size(); ++i)
            if (id_list[i] == pp)
              id_list.erase(id_list.begin() + i);
        }

        // insert entry for new element
        if (_elem_particles.find(elem->unique_id()) != _elem_particles.end())
          _elem_particles[elem->unique_id()].push_back(pp);
        else
          _elem_particles[elem->unique_id()] = {pp};

        // re-assigning the element id
        _md_particles.elem_id[pp] = elem->unique_id();
      }
    }
  }
}

MultiMooseEnum
MDRunBase::getMDQuantities() const
{
  return MultiMooseEnum("vel_x=0 vel_y=1 vel_z=2 force_x=3 force_y=4 force_z=5 charge=6");
}
