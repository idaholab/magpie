/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "NeutronicsSpectrumSamplerFission.h"
#include "MooseMesh.h"

#include <algorithm>

registerMooseObject("MagpieApp", NeutronicsSpectrumSamplerFission);

template <>
InputParameters
validParams<NeutronicsSpectrumSamplerFission>()
{
  InputParameters params = validParams<NeutronicsSpectrumSamplerBase>();
  params.suppressParameter<unsigned int>("L");
  params.set<unsigned int>("L") = 0;
  params.addRequiredCoupledVar("scalar_fluxes", "Scalar fluxes, dimension G.");
  // FIXME: This is not the permanent solution for providing fission cross sections. It must be
  // implemented as material property.
  params.addRequiredParam<std::vector<Real>>(
      "fission_cross_sections", "Fission cross sections. Size = npoints x nisotopes x G.");
  params.addClassDescription("Computes fractional fission rates (Ni * simga_fi * phi) for a "
                             "selection of points.\nUsed for computing PDFs from fission reactions "
                             "to sample PKAs.");
  return params;
}

NeutronicsSpectrumSamplerFission::NeutronicsSpectrumSamplerFission(
    const InputParameters & parameters)
  : NeutronicsSpectrumSamplerBase(parameters)
{
  _nmu = 1;
  _nphi = 1; // can't initialize base class members in initializer

  std::vector<Real> fxs = getParam<std::vector<Real>>("fission_cross_sections");
  if (fxs.size() != _npoints * _I * _G)
    mooseError("fission cross sections must be of length npoints x nisotopes x G");

  // get scalar fluxes
  _scalar_flux.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
    _scalar_flux[g] = &coupledValue("scalar_fluxes", g);

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
NeutronicsSpectrumSamplerFission::totalRecoilRate(unsigned int point_id,
                                                  const std::string & target_isotope) const
{
  Real rate = 0.0;
  auto it = std::find(_target_isotope_names.begin(), _target_isotope_names.end(), target_isotope);
  if (it == _target_isotope_names.end())
    mooseError("Isotope ", target_isotope, "does not exist");
  unsigned int target_isotope_id = std::distance(_target_isotope_names.begin(), it);
  for (unsigned int g = 0; g < _G; ++g)
    rate += _sample_point_data[point_id]({target_isotope_id, g, 0, 0});
  return rate;
}

Real
NeutronicsSpectrumSamplerFission::computeRadiationDamagePDF(unsigned int i,
                                                            unsigned int g,
                                                            unsigned int /*p*/,
                                                            unsigned int /*q*/)
{
  return (*_scalar_flux[g])[_qp] * (*_number_densities[i])[_qp] *
         _fission_cross_section[_current_point][i][g];
}

MultiIndex<Real>
NeutronicsSpectrumSamplerFission::getPDF(unsigned int point_id) const
{
  mooseAssert(_sample_point_data[point_id].size()[2] == 1,
              "RadiationDamageBase: Dimension of mu index is not 1.");
  mooseAssert(_sample_point_data[point_id].size()[3] == 1,
              "RadiationDamageBase: Dimension of phi index is not 1.");
  // the final index of the pdf has dimension 1 so we slice it to return a MI of dimension 2
  return _sample_point_data[point_id].slice(3, 0).slice(2, 0);
}
