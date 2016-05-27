/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMRasterizer.h"
#include "PKAGeneratorBase.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

// custom data load and data store methods for a struct with an std::vector member
template<>
inline void
dataStore(std::ostream & stream, MyTRIMRasterizer::AveragedData & ad, void * context)
{
  dataStore(stream, ad._elements, context);
  dataStore(stream, ad._site_volume, context);
}
template<>
inline void
dataLoad(std::istream & stream, MyTRIMRasterizer::AveragedData & ad, void * context)
{
  dataLoad(stream, ad._elements, context);
  dataLoad(stream, ad._site_volume, context);
}

// custom data load and data store methods for a class with virtual members (vtable pointer must not be (un)serialized)
template<>
inline void
dataStore(std::ostream & stream, MyTRIM_NS::IonBase & pl, void * context)
{
  dataStore(stream, pl._Z, context);
  dataStore(stream, pl._m, context);
  dataStore(stream, pl._E, context);
  dataStore(stream, pl._dir, context);
  dataStore(stream, pl._pos, context);
  dataStore(stream, pl._time, context);
  dataStore(stream, pl.gen, context);
  dataStore(stream, pl.id, context);
  dataStore(stream, pl.tag, context);
  dataStore(stream, pl._Ef, context);
  dataStore(stream, pl.state, context);
}
template<>
inline void
dataLoad(std::istream & stream, MyTRIM_NS::IonBase & pl, void * context)
{
  dataLoad(stream, pl._Z, context);
  dataLoad(stream, pl._m, context);
  dataLoad(stream, pl._E, context);
  dataLoad(stream, pl._dir, context);
  dataLoad(stream, pl._pos, context);
  dataLoad(stream, pl._time, context);
  dataLoad(stream, pl.gen, context);
  dataLoad(stream, pl.id, context);
  dataLoad(stream, pl.tag, context);
  dataLoad(stream, pl._Ef, context);
  dataLoad(stream, pl.state, context);
}

template<>
InputParameters validParams<MyTRIMRasterizer>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Gather the element distribution of the simulation domain for a TRIM binary collision Monte Carlo simulation");
  params.addRequiredCoupledVar("var", "Variables to rasterize");
  params.addCoupledVar("periodic_var", "Optional variables that determines the periodicity. If not supplied the first argument of 'var' will be used.");
  params.addRequiredParam<std::vector<Real> >("M", "Element mass in amu");
  params.addRequiredParam<std::vector<Real> >("Z", "Nuclear charge in e");
  params.addRequiredParam<MaterialPropertyName>("site_volume", "Lattice site volume in nm^3");
  params.addRequiredParam<std::vector<UserObjectName> >("pka_generator", "List of PKA generating user objects");
  MultiMooseEnum setup_options(SetupInterface::getExecuteOptions());
  // we run this object once a timestep
  setup_options = "timestep_begin";
  params.set<MultiMooseEnum>("execute_on") = setup_options;
  return params;
}

MyTRIMRasterizer::MyTRIMRasterizer(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _nvars(coupledComponents("var")),
    _dim(_mesh.dimension()),
    _trim_mass(getParam<std::vector<Real> >("M")),
    _trim_charge(getParam<std::vector<Real> >("Z")),
    _var(_nvars),
    _site_volume_prop(getMaterialProperty<Real>("site_volume")),
    _pka_generator_names(getParam<std::vector<UserObjectName> >("pka_generator")),
    _pka_generators(),
    _periodic(isCoupled("periodic_var") ? coupled("periodic_var", 0) : coupled("var", 0)),
    _last_time(0.0), //TODO: deal with user specified start times!
    _step_end_time(0.0)
{
  for (unsigned int i = 0; i < _nvars; ++i)
    _var[i] = &coupledValue("var", i);

  for (auto && name : _pka_generator_names)
    _pka_generators.push_back(&getUserObjectByName<PKAGeneratorBase>(name));

  if (_nvars == 0)
    mooseError("Must couple variables to MyTRIMRasterier.");

  if (_trim_mass.size() != _nvars)
    mooseError("Parameter 'M' must have as many components as coupled variables.");
  if (_trim_charge.size() != _nvars)
    mooseError("Parameter 'Z' must have as many components as coupled variables.");

  for (unsigned int i = 0; i < _dim; ++i)
  {
    _pbc[i] = _mesh.isRegularOrthogonal() && _mesh.isTranslatedPeriodic(_periodic, i);

    if (_pbc[i])
    {
      _min_dim(i) = _mesh.getMinInDimension(i);
      _max_dim(i) = _mesh.getMaxInDimension(i);
    }
  }
}

bool
MyTRIMRasterizer::executeThisTimestep() const
{
  return true;
}

void
MyTRIMRasterizer::initialize()
{
  _execute_this_timestep = executeThisTimestep();

  // We reset the time of the last run of the BCMC only if the
  // preceeding iteration did converge.
  if (_fe_problem.converged())
    _last_time = _step_end_time;

  if (_execute_this_timestep)
  {
    _material_map.clear();
    _pka_list.clear();

    // Projected time at the end of this step. The total time used to
    // compute the number of PKAs is the time since the end of the last converged
    // timestep in which BCMC ran up to the end of the current timestep.
    _step_end_time = _fe_problem.time() + _fe_problem.dt();
  }
}

void
MyTRIMRasterizer::execute()
{
  // bail out early if not executing this timestep
  if (!_execute_this_timestep)
    return;

  // average element concentrations

  AveragedData average(_nvars);
  Real vol = 0.0;

  // average material data over elements
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    const Real qpvol = _JxW[qp] * _coord[qp];
    vol += qpvol;

    // average compositions on the element
    for (unsigned int i = 0; i < _nvars; ++i)
      average._elements[i] += qpvol * (*_var[i])[qp];

    // average site volume property
    average._site_volume += qpvol * _site_volume_prop[qp];
  }

  // divide by total element volume
  if (vol > 0.0)
  {
    for (unsigned int i = 0; i < _nvars; ++i)
      average._elements[i] /= vol;

    average._site_volume /= vol;
  }

  // store in map
  _material_map[_current_elem->id()] = average;

  // add PKAs for current element
  for (auto && gen : _pka_generators)
    gen->appendPKAs(_pka_list, _step_end_time - _last_time, vol, average);
}

void
MyTRIMRasterizer::threadJoin(const UserObject &y)
{
  // if the map needs to be updated we merge the maps from all threads
  if (_execute_this_timestep)
  {
    const MyTRIMRasterizer & uo = static_cast<const MyTRIMRasterizer &>(y);
    _material_map.insert(uo._material_map.begin(), uo._material_map.end());
    _pka_list.insert(_pka_list.end(), uo._pka_list.begin(), uo._pka_list.end());
  }
}

void
MyTRIMRasterizer::finalize()
{
  // create a one send buffer for use with the libMesh packed range routines
  std::vector<std::string> send_buffers(1);

  // create byte buffers for the streams received from all processors
  std::vector<std::string> recv_buffers;
  recv_buffers.reserve(_app.n_processors());

  // pack the comples datastructures into the string stream
  serialize(send_buffers[0]);

  // broadcast serialized data to and receive from all processors
  _communicator.allgather_packed_range((void *)(nullptr), send_buffers.begin(), send_buffers.end(),
                                       std::back_inserter(recv_buffers));

  // unpack the received data and merge it into the local data structures
  deserialize(recv_buffers);
}

const std::vector<Real> &
MyTRIMRasterizer::material(const Elem * elem) const
{
  auto i = _material_map.find(elem->id());

  // there should be data for every element in the mesh
  if (i == _material_map.end())
    mooseError("Element not found in material map.");

  return i->second._elements;
}

Real
MyTRIMRasterizer::siteVolume(const Elem * elem) const
{
  auto i = _material_map.find(elem->id());

  // there should be data for every element in the mesh
  if (i == _material_map.end())
    mooseError("Element not found in material map.");

  return i->second._site_volume;
}

Point
MyTRIMRasterizer::periodicPoint(const Point & pos) const
{
  // point to sample the material at
  Point p(pos(0), pos(1), _dim == 2 ? 0.0 : pos(2));

  // apply periodic boundary conditions
  for (unsigned int i = 0; i < _dim; ++i)
    if (_pbc[i])
    {
      const Real width = _max_dim(i) - _min_dim(i);
      p(i) -= std::floor((p(i) - _min_dim(i)) / width) * width;
    }

  return p;
}

void
MyTRIMRasterizer::serialize(std::string & serialized_buffer)
{
  // stream for serializing the _material_map and _pka_list data structure to a byte stream
  std::ostringstream oss;
  dataStore(oss, _material_map, this);
  dataStore(oss, _pka_list, this);

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
}

void
MyTRIMRasterizer::deserialize(std::vector<std::string> & serialized_buffers)
{
  mooseAssert(serialized_buffers.size() == _app.n_processors(), "Unexpected size of serialized_buffers: " << serialized_buffers.size());

  // The input string stream used for deserialization
  std::istringstream iss;

  // Loop over all datastructures for all procerssors to perfrom the gather operation
  for (unsigned int rank = 0; rank < serialized_buffers.size(); ++rank)
  {
    // skip the current processor (its data is already in the strutures)
    if (rank == processor_id())
      continue;

    // populate the stream with a new buffer and reset stream state
    iss.str(serialized_buffers[rank]);
    iss.clear();

    // Load the communicated data into temporary structures
    MaterialMap other_material_map;
    dataLoad(iss, other_material_map, this);
    std::vector<MyTRIM_NS::IonBase> other_pka_list;
    dataLoad(iss, other_pka_list, this);

    // merge the data in with the current processor's data
    _material_map.insert(other_material_map.begin(), other_material_map.end());

    // We are not yet merging the PKA lists but let each processor work on its own list.
    // List combining will be used in the furture to enable better load balancing.
    // _pka_list.insert(other_pka_list.begin(), other_pka_list.end());
  }
}
