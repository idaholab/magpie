/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "IsotopeRecoilRate.h"
#include "MyTRIMRasterizer.h"
#include "NeutronicsSpectrumSamplerBase.h"
#include "mytrim/ion.h"

registerMooseObject("MagpieApp", IsotopeRecoilRate);

InputParameters
IsotopeRecoilRate::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<std::string>(
      "target_isotope", "The isotope name that you want to get the total recoil rate for");
  params.addRequiredParam<unsigned int>("point_id", "The index of the point in neutronics_sampler");
  params.addRequiredParam<UserObjectName>(
      "neutronics_sampler", "The neutronics sampler object that the data is retrieved from");
  params.addParam<PostprocessorName>(
      "scaling_factor", 1, "A scaling factor multiplying the isotope recoil rate");
  params.addClassDescription("Gets the total recoil rate from target_isotope at point point_id "
                             "contained in the neutronics_sampler");
  return params;
}

IsotopeRecoilRate::IsotopeRecoilRate(const InputParameters & params)
  : GeneralPostprocessor(params),
    _target_isotope(getParam<std::string>("target_isotope")),
    _point_id(getParam<unsigned int>("point_id")),
    _neutronics_sampler(getUserObject<NeutronicsSpectrumSamplerBase>("neutronics_sampler")),
    _scaling_factor(getPostprocessorValue("scaling_factor"))
{
  if (_neutronics_sampler.getNumberOfPoints() < _point_id)
    mooseError("The provided neutronics sampler object only has",
               _neutronics_sampler.getNumberOfPoints(),
               " points");

  if (!_neutronics_sampler.hasIsotope(_target_isotope))
    mooseError("Target isotope ", _target_isotope, " not preset in neutronics sampler object");
}

Real
IsotopeRecoilRate::getValue() const
{
  return _scaling_factor * _neutronics_sampler.totalRecoilRate(_point_id, _target_isotope);
}
