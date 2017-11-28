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

#include "IsotopeRecoilRateSampler.h"
#include "NeutronicsSpectrumSamplerBase.h"

// MOOSE includes
#include "MooseMesh.h"
#include "MooseVariable.h"

#include "libmesh/mesh_tools.h"

template <>
InputParameters
validParams<IsotopeRecoilRateSampler>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addRequiredParam<std::string>("target_isotope", "The isotope name that you want to get the total recoil rate for");
  params.addRequiredParam<std::vector<unsigned int>>("point_ids", "The indices of the points in neutronics_sampler");
  params.addRequiredParam<UserObjectName>("neutronics_sampler", "The neutronics sampler object that the data is retrieved from");
  params.addClassDescription("Gets the total recoil rate from target_isotope at points provided in point_id contained in the neutronics_sampler");
  return params;
}

IsotopeRecoilRateSampler::IsotopeRecoilRateSampler(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
  _target_isotope(getParam<std::string>("target_isotope")),
  _point_ids(getParam<std::vector<unsigned int>>("point_ids")),
  _neutronics_sampler(getUserObject<NeutronicsSpectrumSamplerBase>("neutronics_sampler")),
  _recoil_rates(declareVector("recoil_rates"))
{
  _recoil_rates.assign(_point_ids.size(), 0);
  for (auto & p : _point_ids)
    if (_neutronics_sampler.getNumberOfPoints() < p)
      mooseError("The provided neutronics sampler object only has", _neutronics_sampler.getNumberOfPoints(), " points but point id ", p, " is requested");

  if (!_neutronics_sampler.hasIsotope(_target_isotope))
    mooseError("Target isotope ", _target_isotope, " not preset in neutronics sampler object");
}

void
IsotopeRecoilRateSampler::execute()
{
  for (unsigned int j = 0; j < _point_ids.size(); ++j)
    _recoil_rates[j] = _neutronics_sampler.totalRecoilRate(_point_ids[j], _target_isotope);
}
