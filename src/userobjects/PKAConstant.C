/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAConstant.h"

registerMooseObject("MagpieApp", PKAConstant);

InputParameters
PKAConstant::validParams()
{
  InputParameters params = PKAEmpiricalBase::validParams();
  params.addClassDescription("PKAs with constant mass, charge, energy, and rate");
  params.addParam<Real>(
      "pka_rate",
      1e-8,
      "PKA rate per unit volume (uses mesh units defined in the rasterizer and moose time units)");
  params.addRequiredParam<Real>("Z", "PKA nuclear charge");
  params.addRequiredParam<Real>("m", "PKA mass in amu");
  params.addRequiredParam<Real>("E", "PKA energy in eV");
  return params;
}

PKAConstant::PKAConstant(const InputParameters & parameters)
  : PKAEmpiricalBase(parameters),
    _pka_rate(getParam<Real>("pka_rate")),
    _Z(getParam<Real>("Z")),
    _m(getParam<Real>("m")),
    _E(getParam<Real>("E"))
{
}
