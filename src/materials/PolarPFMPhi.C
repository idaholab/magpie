/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPFMPhi.h"

registerMooseObject("MagpieApp", PolarPFMPhi);

InputParameters
PolarPFMPhi::validParams()
{
  InputParameters params = DerivativeParsedMaterialHelper::validParams();
  params.addClassDescription("Material property for phi with a beta21 prefactor");
  params.addRequiredCoupledVar("upsilon", "Upsilon order parameter");
  params.addRequiredParam<Real>("a_phi", "Interpolation coefficient a_phi");
  params.addRequiredParam<Real>("a0", "Interpolation coefficient a0");
  params.addRequiredParam<Real>("beta21",
                                "Gradient energy coefficient between solid 2 and solid 1");
  return params;
}

PolarPFMPhi::PolarPFMPhi(const InputParameters & parameters)
  : DerivativeParsedMaterialHelper(parameters),
    _upsilon("upsilon"),
    _a_phi(getParam<Real>("a_phi")),
    _a0(getParam<Real>("a0")),
    _beta21(getParam<Real>("beta21"))
{
  // interpolating function (9) from Physical Review B 89, 184102 (2014)
  EBTerm phi = _a_phi * _upsilon * _upsilon -
               2 * (_a_phi - 2 * (1 - _a0)) * _upsilon * _upsilon * _upsilon +
               (_a_phi - 3 * (1 - _a0)) * _upsilon * _upsilon * _upsilon * _upsilon + _a0;

  functionParse(_beta21 * phi);
}
