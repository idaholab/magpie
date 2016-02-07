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

template<>
InputParameters validParams<MyTRIMRasterizer>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Gather the element distribution of the simulation domain for a TRIM binary collision Monte Carlo simulation");
  params.addCoupledVar("var", "Variables to rasterize");
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
    _pka_generators(_pka_generator_names.size()),
    _periodic(coupled("var", 0)),
    _last_time(0.0), //TODO: deal with user specified start times!
    _step_end_time(0.0)
{
  for (unsigned int i = 0; i < _nvars; ++i)
    _var[i] = &coupledValue("var", i);

  for (unsigned int i = 0; i < _pka_generator_names.size(); ++i)
    _pka_generators[i] = &getUserObjectByName<PKAGeneratorBase>(_pka_generator_names[i]);

  if (_nvars == 0)
    mooseError("Must couple variables to MyTRIMRasterier.");

  if (_trim_mass.size() != _nvars)
    mooseError("Parameter 'M' must have as many components as coupled variables.");
  if (_trim_charge.size() != _nvars)
    mooseError("Parameter 'Z' must have as many components as coupled variables.");

  if (_app.n_processors() > 1)
    mooseError("Parallel communication is not yet implemented.");

  for (unsigned int i = 0; i < _dim; ++i)
  {
    _pbc[i] = _mesh.isTranslatedPeriodic(_periodic, i);

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
  for (unsigned int i = 0; i < _pka_generators.size(); ++i)
    _pka_generators[i]->appendPKAs(_pka_list, _step_end_time - _last_time, vol, average);
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
  // _communicator.set_union(_material_map);
}

const std::vector<Real> &
MyTRIMRasterizer::material(const Elem * elem) const
{
  MaterialMap::const_iterator i = _material_map.find(elem->id());

  // there should be data for every element in the mesh
  if (i == _material_map.end())
    mooseError("Element not found in material map.");

  return i->second._elements;
}

Real
MyTRIMRasterizer::siteVolume(const Elem * elem) const
{
  MaterialMap::const_iterator i = _material_map.find(elem->id());

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
