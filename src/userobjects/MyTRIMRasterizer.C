/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMRasterizer.h"

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
  MultiMooseEnum setup_options(SetupInterface::getExecuteOptions());
  // we run this object once a timestep
  setup_options = "timestep_begin";
  params.set<MultiMooseEnum>("execute_on") = setup_options;
  return params;
}

MyTRIMRasterizer::MyTRIMRasterizer(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _nvars(coupledComponents("var")),
    _trim_mass(getParam<std::vector<Real> >("M")),
    _trim_charge(getParam<std::vector<Real> >("Z")),
    _var(_nvars)
{
  for (unsigned int i = 0; i < _nvars; ++i)
    _var[i] = &coupledValue("var", i);

  if (_nvars == 0)
    mooseError("Must couple variables to MyTRIMRasterier.");

  if (_trim_mass.size() != _nvars)
    mooseError("Parameter 'M' must have as many components as coupled variables.");
  if (_trim_charge.size() != _nvars)
    mooseError("Parameter 'Z' must have as many components as coupled variables.");

  _periodic = coupled("var", 0);

  if (_app.n_processors() > 1)
    mooseError("Parallel communication is not yet implemented. Waiting on libmesh/#748.");
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

  if (_execute_this_timestep)
    _material_map.clear();
}

void
MyTRIMRasterizer::execute()
{
  // bail out early if not executing this timestep
  if (!_execute_this_timestep)
    return;

  // average element concentrations
  std::vector<Real> elements(_nvars, 0.0);
  Real vol = 0.0;

  // average material data over elements
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    const Real qpvol = _JxW[qp] * _coord[qp];
    vol += qpvol;
    for (unsigned int i = 0; i < _nvars; ++i)
      elements[i] += qpvol * (*_var[i])[qp];
  }

  // divide by total element volume
  if (vol > 0.0)
    for (unsigned int i = 0; i < _nvars; ++i)
      elements[i] /= vol;

  // store in map
  _material_map[_current_elem->id()] = elements;
}

void
MyTRIMRasterizer::threadJoin(const UserObject &y)
{
  // if the map needs to be updated we merge the maps from all threads
  if (_execute_this_timestep)
  {
    const MyTRIMRasterizer & uo = static_cast<const MyTRIMRasterizer &>(y);
    _material_map.insert(uo._material_map.begin(), uo._material_map.end());
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

  return i->second;
}
