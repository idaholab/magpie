/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPFMGradient.h"

registerMooseObject("MagpieApp", PolarPFMGradient);

InputParameters
PolarPFMGradient::validParams()
{
  InputParameters params = KernelValue::validParams();
  params.addClassDescription("Gradient energy term in the polar phase field model");
  params.addRequiredParam<MaterialPropertyName>(
      "prop", "Material property of which the kernel variable derivative will be taken of");
  params.addCoupledVar("v", "Coupled order parameter");
  return params;
}

PolarPFMGradient::PolarPFMGradient(const InputParameters & parameters)
  : DerivativeMaterialInterface<KernelValue>(parameters),
    _grad_v(coupledGradient("v")),
    _v_name(getVar("v", 0)->name()),
    _v_var(coupled("v")),
    _dpropdu(getMaterialPropertyDerivative<Real>("prop", _var.name())),
    _d2propdu2(getMaterialPropertyDerivative<Real>("prop", _var.name(), _var.name())),
    _d2propdudv(getMaterialPropertyDerivative<Real>("prop", _var.name(), _v_name))
{
}

Real
PolarPFMGradient::precomputeQpResidual()
{
  return -0.5 * _dpropdu[_qp] * _grad_v[_qp].norm_sq();
}

Real
PolarPFMGradient::precomputeQpJacobian()
{
  return -0.5 * _d2propdu2[_qp] * _phi[_j][_qp] * _grad_v[_qp].norm_sq();
}

Real
PolarPFMGradient::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_var)
    return -0.5 *
           (_d2propdudv[_qp] * _phi[_j][_qp] * _grad_v[_qp].norm_sq() +
            _dpropdu[_qp] * 2.0 * _grad_v[_qp] * _grad_phi[_j][_qp]) *
           _test[_i][_qp];

  return 0.0;
}
