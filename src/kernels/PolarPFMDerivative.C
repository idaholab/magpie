/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPFMDerivative.h"

registerMooseObject("MagpieApp", PolarPFMDerivative);

InputParameters
PolarPFMDerivative::validParams()
{
  InputParameters params = KernelValue::validParams();
  params.addClassDescription("Bulk energy derivative term in the polar phase field model");
  params.addRequiredParam<MaterialPropertyName>(
      "F", "Material property of which the kernel variable derivative will be taken of");
  return params;
}

PolarPFMDerivative::PolarPFMDerivative(const InputParameters & parameters)
  : DerivativeMaterialInterface<KernelValue>(parameters),
    _dpropdu(getMaterialPropertyDerivative<Real>("F", _var.name())),
    _d2propdu2(getMaterialPropertyDerivative<Real>("F", _var.name(), _var.name()))

{
}

Real
PolarPFMDerivative::precomputeQpResidual()
{
  return _dpropdu[_qp];
}

Real
PolarPFMDerivative::precomputeQpJacobian()
{
  return _d2propdu2[_qp] * _phi[_j][_qp];
}
