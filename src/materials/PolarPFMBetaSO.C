/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPFMBetaSO.h"

registerMooseObject("MagpieApp", PolarPFMBetaSO);

InputParameters
PolarPFMBetaSO::validParams()
{
  InputParameters params = DerivativeParsedMaterialHelper::validParams();
  params.addClassDescription("Material property betaSO");
  params.addRequiredCoupledVar("theta", "Theta order parameter");
  params.addRequiredParam<Real>("a_beta", "Interpolation coefficient a_beta");
  params.addRequiredParam<Real>("beta10", "Gradient energy coefficient between solid 1 and melt");
  params.addRequiredParam<Real>("beta20", "Gradient energy coefficient between solid 2 and melt");
  return params;
}

PolarPFMBetaSO::PolarPFMBetaSO(const InputParameters & parameters)
  : DerivativeParsedMaterialHelper(parameters),
    _theta("theta"),
    _a_beta(getParam<Real>("a_beta")),
    _beta10(getParam<Real>("beta10")),
    _beta20(getParam<Real>("beta20"))
{
  // All equation numbers from Physical Review B 89, 184102 (2014)

  // interpolating function (8)
  EBTerm q = _a_beta * _theta * _theta - 2 * (_a_beta - 2) * _theta * _theta * _theta +
             (_a_beta - 3) * _theta * _theta * _theta * _theta;

  // solid-melt gradient energy coefficient (7)
  EBTerm betaSO = _beta10 + (_beta20 - _beta10) * q;

  functionParse(betaSO);
}
