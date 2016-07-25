/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifdef RATTLESNAKE_ENABLED
#include "RadiationDamageBase.h"
#include "MooseMesh.h"
#include "YakxsUtilities.h"

// C++ includes
#include <sstream>
#include <algorithm>
#include <limits>

template<>
InputParameters validParams<RadiationDamageBase>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredParam<std::vector<std::string> >("target_isotope_names", "The list of target isotope names e.g. U235.");
  params.addRequiredCoupledVar("number_densities", "Number densities for each isotope.");
  params.addRequiredParam<std::vector<Real> >("energy_group_boundaries", "The energy group boundaries ommitting E = 0.0. Units are MeV.");
  params.addRequiredParam<unsigned int>("L", "The order up to which angular moments of the PKA distribution are computed.");
  params.addRequiredParam<std::vector<Point> >("points", "The points where you want to evaluate the variables");
  params.addClassDescription("Primary Knock-on Atom (PKA) user object base class. Computes PKA distributions at a selection of points.");
  return params;
}

RadiationDamageBase::RadiationDamageBase(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _target_isotope_names(getParam<std::vector<std::string> >("target_isotope_names")),
    _energy_group_boundaries(getParam<std::vector<Real> >("energy_group_boundaries")),
    _I(_target_isotope_names.size()),
    _G(_energy_group_boundaries.size()),
    _L(getParam<unsigned int>("L")),
    _points(getParam<std::vector<Point> >("points")),
    _npoints(_points.size()),
    _qp_is_cached(false)
{
  // check input dimensions
  if (coupledComponents("number_densities") != _I)
    mooseError("ZAID and number_densities must have the same length.");

  // get the number density variables
  _number_densities.resize(_I);
  for (unsigned int i = 0; i < _I; ++i)
    _number_densities[i] = & coupledValue("number_densities", i);

  // check if isotope names are valid, NOTE: error handling is delegated to Yakxs::Utilities
  unsigned int Z, A;
  for (unsigned int i = 0; i < _I; ++i)
    YAKXS::Utility::getAZFromIsotopeName(_target_isotope_names[i], A, Z);
  _point_element.resize(_npoints);
  _qp_cache.resize(_npoints);
}

void
RadiationDamageBase::execute()
{
  if (std::find(_point_element.begin(), _point_element.end(), _current_elem->id()) != _point_element.end())
  {
    for (_current_point = 0; _current_point < _npoints; _current_point++)
    {
      if (_point_element[_current_point] != _current_elem->id()) continue;

      // NOTE: PKA is evaluated at a point that is not necessarily a qp but material
      // props only live on qps. This is a cheap and dirty solution evaluating the
      // PKA at the closest quadrature point.
      if (!_qp_is_cached)
      {
        Real min_dsq = (_q_point[0] - _points[_current_point]).norm_sq();
        unsigned int min_qp = 0;
        for (unsigned int qp = 1; qp < _q_point.size(); qp++)
        {
          Real dsq = (_q_point[qp] - _points[_current_point]).norm_sq();
          if (dsq < min_dsq)
          {
            min_dsq = dsq;
            min_qp = qp;
          }
        }
        _qp_cache[_current_point] = min_qp;
      }

      // Set current qp to cached min_qp
      _qp = _qp_cache[_current_point];

      // call back before computing PKA distribution for caching things that
      // don't change
      preComputePKA();

      // Evalute the PKA distribution at min_qp
      for (unsigned int i = 0; i < _I; ++i)
        for (unsigned int g = 0; g < _G; ++g)
          for (unsigned int p = 0; p < _nSH; ++p)
            _sample_point_data[_current_point][i][g][p] = computePKA(i, g, p);
    }
  }
}

void
RadiationDamageBase::preComputePKA()
{
}

void
RadiationDamageBase::meshChanged()
{
  // Rebuild point_locator & find the right element for each point
  UniquePtr<PointLocatorBase> point_locator = PointLocatorBase::build(TREE_ELEMENTS, _mesh.getMesh());
  for (unsigned int j = 0; j < _npoints; j++)
    _point_element[j] = (*point_locator)(_points[j])->id();
  _qp_is_cached = false;
}

void
RadiationDamageBase::initialize()
{
  // allocate PKA distribution
  // NOTE: Needs to be delayed to initialize because _nSH is set in derived class
  _sample_point_data.resize(_npoints);
  for (unsigned j = 0; j < _npoints; ++j)
  {
    _sample_point_data[j].resize(_I);
    for (unsigned int i = 0; i < _I; ++i)
    {
      _sample_point_data[j][i].resize(_G);
      for (unsigned int g = 0; g < _G; ++g)
        _sample_point_data[j][i][g].resize(_nSH);
    }
  }

  meshChanged();
}

void
RadiationDamageBase::finalize()
{
  for (unsigned j = 0; j < _npoints; ++j)
    for (unsigned int i = 0; i < _I; ++i)
      for (unsigned int g = 0; g < _G; ++g)
        for (unsigned int p = 0; p < _nSH; ++p)
          gatherSum(_sample_point_data[j][i][g][p]);

  // set _qp_is_cached flag to true. Actually we only need to do this if
  // it was false but this way we can save an if statement
  _qp_is_cached = true;
}

void
RadiationDamageBase::threadJoin(const UserObject & y)
{
  const RadiationDamageBase & uo = static_cast<const RadiationDamageBase &>(y);
  for (unsigned j = 0; j < _npoints; ++j)
    for (unsigned int i = 0; i < _I; ++i)
      for (unsigned int g = 0; g < _G; ++g)
        for (unsigned int p = 0; p < _nSH; ++p)
          _sample_point_data[j][i][g][p] += uo._sample_point_data[j][i][g][p];
}

MultiIndex<Real>
getPDF(unsigned int point_id) const
{
  return _sample_point_data[point_id];
}

Real
getMagnitude(unsigned int point_id) const
{
  return _magnitude[point_id];
}

std::vector<unsigned int>
getZAIDs(unsigned int point_id) const
{
  return _zaids[point_id];
}

std::vector<Real>
getEnergies(unsigned int point_id) const
{
  return _energies[point_id];
}

#endif //RATTLESNAKE_ENABLED
