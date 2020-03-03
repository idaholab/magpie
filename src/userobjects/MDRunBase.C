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
  dataStore(stream, pl.id, context);
  dataStore(stream, pl.elem_id, context);
  dataStore(stream, pl.properties, context);
}

template <>
inline void
dataLoad(std::istream & stream, MDRunBase::MDParticles & pl, void * context)
{
  dataLoad(stream, pl.pos, context);
  dataLoad(stream, pl.id, context);
  dataLoad(stream, pl.elem_id, context);
  dataLoad(stream, pl.properties, context);
}

InputParameters
MDRunBase::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  params.suppressParameter<ExecFlagEnum>("execute_on");
  params.addParam<MultiMooseEnum>("md_particle_properties",
                                  MDRunBase::mdParticleProperties(),
                                  "Properties of MD particles to be obtained and stored.");
  params.addRangeCheckedParam<Real>(
      "max_granular_radius",
      0,
      "max_granular_radius>=0",
      "Maximum radius of granular particles. This can be important for partitioning.");
  params.addClassDescription(
      "Base class for execution of coupled molecular dynamics MOOSE calculations.");
  return params;
}

MDRunBase::MDRunBase(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _properties(getParam<MultiMooseEnum>("md_particle_properties")),
    _granular(_properties.contains("radius")),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension()),
    _nproc(_app.n_processors())
{
  _max_granular_radius = getParam<Real>("max_granular_radius");
  // if the calculation is granular, max_particle_radius must be set
  if (_granular && !parameters.isParamSetByUser("max_granular_radius"))
    paramError("max_granular_radius", "max_granular_radius must be set for granular calculations");

  // set up the map from property ID to index
  _md_particles._prop_size = _properties.size();
  for (unsigned int j = 0; j < _properties.size(); ++j)
    _md_particles._map_props[_properties.get(j)] = j;

  // _md_particles._r_index is a short-cut for the "radius" index
  if (_granular)
    _md_particles._r_index = propIndex("radius");
}

void
MDRunBase::initialSetup()
{
  _bbox = MeshTools::create_processor_bounding_box(_mesh, processor_id());

  // inflate bounding box
  for (unsigned int d = 0; d < _dim; ++d)
  {
    _bbox.first(d) -= _max_granular_radius;
    _bbox.second(d) += _max_granular_radius;
  }
}

void
MDRunBase::timestepSetup()
{
  // update/init the processor bounding box
  _bbox = MeshTools::create_processor_bounding_box(_mesh, processor_id());

  // inflate bounding box
  for (unsigned int d = 0; d < _dim; ++d)
  {
    _bbox.first(d) -= _max_granular_radius;
    _bbox.second(d) += _max_granular_radius;
  }

  // callback for updating md particle list
  updateParticleList();
}

void
MDRunBase::updateKDTree()
{
  _kd_tree = libmesh_make_unique<KDTree>(_md_particles.pos, 50);
}

void
MDRunBase::elemParticles(unique_id_type elem_id, std::vector<unsigned int> & elem_particles) const
{
  if (_elem_particles.find(elem_id) != _elem_particles.end())
    elem_particles = _elem_particles.find(elem_id)->second;
  else
    elem_particles = {};
}

void
MDRunBase::granularElementVolumes(unique_id_type elem_id,
                                  std::vector<std::pair<unsigned int, Real>> & gran_vol) const
{
  mooseAssert(_granular,
              "Radius must be provided as MD property to allow granular volume computation.");
  if (_elem_granular_volumes.find(elem_id) != _elem_granular_volumes.end())
    gran_vol = _elem_granular_volumes.find(elem_id)->second;
  else
    gran_vol = {};
}

Real
MDRunBase::particleProperty(unsigned int j, unsigned int prop_id) const
{
  // ensure that entry exists
  mooseAssert(j < _md_particles.properties.size(),
              "Particle index " << j << " not found in _md_particles. properties vector has length "
                                << _md_particles.properties.size());
  return _md_particles.properties[j][propIndex(prop_id)];
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

void
MDRunBase::updateElementGranularVolumes()
{
  // clear the granular volume map
  _elem_granular_volumes.clear();

  /*
     // Ideally we want to compute _max_granular_radius automatically
     // but the value is needed for setting bounding boxes to partition
     // the particles on the processors. For now, it's an input parameter
     // until a path forward is determined.
     // The commented lines compute the largest granular radius
  _max_granular_radius = 0;
  for (auto & p : _md_particles.properties)
    if (p[7] > _max_granular_radius)
      _max_granular_radius = p[7];
  */
  /// do a sanity check _max_granular_radius
  Real mgr = 0;
  for (auto & p : _md_particles.properties)
    if (p[_md_particles._r_index] > mgr)
      mgr = p[_md_particles._r_index];
  if (mgr > _max_granular_radius)
    mooseError("Granular particle with radius: ",
               mgr,
               " exceeds max_granular_radius: ",
               _max_granular_radius);

  /// loop over all local elements
  ConstElemRange * active_local_elems = _mesh.getActiveLocalElementRange();
  for (const auto & elem : *active_local_elems)
  {
    // find all points within an inflated bounding box
    std::vector<std::pair<std::size_t, Real>> indices_dist;
    BoundingBox bbox = elem->loose_bounding_box();
    Point center = 0.5 * (bbox.min() + bbox.max());

    // inflate the search sphere by the maximum granular radius
    Real radius = (bbox.max() - center).norm() + _max_granular_radius;
    _kd_tree->radiusSearch(center, radius, indices_dist);

    // prepare _elem_granular_candidates entry
    _elem_granular_volumes[elem->unique_id()] = {};

    // construct this element's overlap object
    ElemType t = elem->type();
    OVERLAP::Hexahedron hex = overlapUnitHex();
    OVERLAP::Tetrahedron tet = overlapUnitTet();
    if (t == HEX8)
      hex = overlapHex(elem);
    else if (t == TET4)
      tet = overlapTet(elem);
    else
      mooseError("Element type ", t, "not implemented");

    // loop through all MD particles that the search turned up and test overlap
    for (unsigned int j = 0; j < indices_dist.size(); ++j)
    {
      // construct OVERLAP::sphere object from MD granular particle
      unsigned int k = indices_dist[j].first;
      OVERLAP::Sphere sph(OVERLAP::vector_t{_md_particles.pos[k](0),
                                            _md_particles.pos[k](1),
                                            _md_particles.pos[k](2)},
                          _md_particles.properties[k][_md_particles._r_index]);

      // compute the overlap
      Real ovlp = 0.0;
      if (t == HEX8)
        ovlp = OVERLAP::overlap(sph, hex);
      else if (t == TET4)
        ovlp = OVERLAP::overlap(sph, tet);

      // if the overlap is larger than 0, make entry in _elem_granular_volumes
      if (ovlp > 0.0)
        _elem_granular_volumes[elem->unique_id()].push_back(std::pair<unsigned int, Real>(k, ovlp));
    }
  }
}

MultiMooseEnum
MDRunBase::mdParticleProperties()
{
  return MultiMooseEnum("vel_x=0 vel_y=1 vel_z=2 force_x=3 force_y=4 force_z=5 charge=6 radius=7");
}

unsigned int
MDRunBase::propIndex(unsigned int prop_id) const
{
  auto it = _md_particles._map_props.find(prop_id);
  if (it == _md_particles._map_props.end())
    mooseError("Property id ", prop_id, " is not present in _map_props map.");
  return it->second;
}

unsigned int
MDRunBase::propIndex(const std::string & prop_name) const
{
  unsigned int prop_id = mdParticleProperties().find(prop_name)->id();
  return propIndex(prop_id);
}

OVERLAP::Hexahedron
MDRunBase::overlapHex(const Elem * elem) const
{
  Point p;
  p = elem->point(0);
  OVERLAP::vector_t v0{p(0), p(1), p(2)};
  p = elem->point(1);
  OVERLAP::vector_t v1{p(0), p(1), p(2)};
  p = elem->point(2);
  OVERLAP::vector_t v2{p(0), p(1), p(2)};
  p = elem->point(3);
  OVERLAP::vector_t v3{p(0), p(1), p(2)};
  p = elem->point(4);
  OVERLAP::vector_t v4{p(0), p(1), p(2)};
  p = elem->point(5);
  OVERLAP::vector_t v5{p(0), p(1), p(2)};
  p = elem->point(6);
  OVERLAP::vector_t v6{p(0), p(1), p(2)};
  p = elem->point(7);
  OVERLAP::vector_t v7{p(0), p(1), p(2)};
  return OVERLAP::Hexahedron{v0, v1, v2, v3, v4, v5, v6, v7};
}

OVERLAP::Hexahedron
MDRunBase::overlapUnitHex() const
{
  OVERLAP::vector_t v0{-1, -1, -1};
  OVERLAP::vector_t v1{1, -1, -1};
  OVERLAP::vector_t v2{1, 1, -1};
  OVERLAP::vector_t v3{-1, 1, -1};
  OVERLAP::vector_t v4{-1, -1, 1};
  OVERLAP::vector_t v5{1, -1, 1};
  OVERLAP::vector_t v6{1, 1, 1};
  OVERLAP::vector_t v7{-1, 1, 1};
  return OVERLAP::Hexahedron{v0, v1, v2, v3, v4, v5, v6, v7};
}

OVERLAP::Tetrahedron
MDRunBase::overlapTet(const Elem * elem) const
{
  Point p;
  p = elem->point(0);
  OVERLAP::vector_t v0{p(0), p(1), p(2)};
  p = elem->point(1);
  OVERLAP::vector_t v1{p(0), p(1), p(2)};
  p = elem->point(2);
  OVERLAP::vector_t v2{p(0), p(1), p(2)};
  p = elem->point(3);
  OVERLAP::vector_t v3{p(0), p(1), p(2)};
  return OVERLAP::Tetrahedron{v0, v1, v2, v3};
}

OVERLAP::Tetrahedron
MDRunBase::overlapUnitTet() const
{
  OVERLAP::vector_t v0{0, 0, 0};
  OVERLAP::vector_t v1{1, 0, 0};
  OVERLAP::vector_t v2{0, 1, 0};
  OVERLAP::vector_t v3{0, 0, 1};
  return OVERLAP::Tetrahedron{v0, v1, v2, v3};
}
