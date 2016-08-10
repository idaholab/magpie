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
#include "NeutronicsSpectrumSamplerFission.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<NeutronicsSpectrumSamplerFission>()
{
  InputParameters params = validParams<NeutronicsSpectrumSamplerBase>();
  params.suppressParameter<unsigned int>("L");
  params.set<unsigned int>("L") = 0;
  params.addRequiredCoupledVar("scalar_fluxes", "Scalar fluxes, dimension G.");
  // FIXME: This is not the permanent solution for providing fission cross sections. It must be implemented as material
  // property.
  params.addRequiredParam<std::vector<Real> >("fission_cross_sections", "Fission cross sections. Size = npoints x nisotopes x G.");
  params.addClassDescription("Computes fractional fission rates (Ni * simga_fi * phi) for a selection of points.\nUsed for computing PDFs from fission reactions to sample PKAs.");
  return params;
}

NeutronicsSpectrumSamplerFission::NeutronicsSpectrumSamplerFission(const InputParameters & parameters) :
    NeutronicsSpectrumSamplerBase(parameters)
{
  _nSH = 1; // can't initialize base class members in initializer
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
NeutronicsSpectrumSamplerFission::computeRadiationDamagePDF(unsigned int i, unsigned int g, unsigned int /*p*/)
{
  return (*_scalar_flux[g])[_qp] * (*_number_densities[i])[_qp] * _fission_cross_section[_current_point][i][g];
}

MultiIndex<Real>
NeutronicsSpectrumSamplerFission::getPDF(unsigned int point_id) const
{
  mooseAssert(_sample_point_data[point_id].size()[2] == 1, "RadiationDamageBase: Dimension of last index is not 1.");
  // the final index of the pdf has dimension 1 so we slice it to return a MI of dimension 2
  return _sample_point_data[point_id].slice(2, 0);
}
