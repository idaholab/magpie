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
#include "RadiationDamageFission.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<RadiationDamageFission>()
{
  InputParameters params = validParams<RadiationDamageBase>();
  params.suppressParameter<unsigned int>("L");
  params.set<unsigned int>("L") = 0;
  params.addRequiredCoupledVar("scalar_fluxes", "Scalar fluxes, dimension G.");
  // FIXME: This is not the permanent solution for providing fission cross sections. It must be implemented as material
  // property.
  params.addRequiredParam<std::vector<Real> >("fission_cross_sections", "Fission cross sections. Size = npoints x nisotopes x G.");
  params.addClassDescription("Computes fractional fission rates (Ni * simga_fi * phi) for a selection of points.\nUsed for computing PDFs from fission reactions to sample PKAs.");
  return params;
}

RadiationDamageFission::RadiationDamageFission(const InputParameters & parameters) :
    RadiationDamageBase(parameters),
    _nSH = 1
{
  std::vector<Real> fxs = getParam<std::vector<Real> >("fission_cross_sections");
  if (fxs.size() != _npoints * _I * _G)
    mooseError("fission cross sections must be of length npoints x nisotopes x G");

  // get scalar fluxes
  _scalar_flux.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
    _scalar_flux[g] = & coupledValue("scalar_fluxes", g);

  // allocate and assign fission cross sections
  unsigned int p = 0;
  _fission_cross_section.resize(_npoints);
  for (unsigned j = 0; j < _npoints; ++j)
  {
    _fission_cross_section[j].resize(_I);
    for (unsigned int i = 0; i < _I; ++i)
    {
      _fission_cross_section[j][i].resize(_G);
      for (unsigned int g = 0; g < _G; ++g)
        _fission_cross_section[j][i][g] = fxs[p++];
    }
  }
}

Real
RadiationDamageFission::computeRadiationDamagePDF(unsigned int i, unsigned int g, unsigned int /*p*/)
{
  return (*_scalar_flux[g])[_qp] * (*_number_densities[i])[_qp] * _fission_cross_section[_current_point][i][g];
}

#endif //RATTLESNAKE_ENABLED
