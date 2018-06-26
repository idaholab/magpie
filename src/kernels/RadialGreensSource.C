/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "RadialGreensSource.h"

registerMooseObject("MagpieApp", RadialGreensSource);

template <>
InputParameters
validParams<RadialGreensSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Apply the convolution from a RadialGreensConvolution object to a non-linear variable");
  params.addRequiredParam<UserObjectName>("convolution", "RadialGreensConvolution user object");
  params.addCoupledVar("c",
                       "coupled concentration variable. If this is not supplied the kernel "
                       "variable is used. This is useful in conjunction with "
                       "CoupledTimeDerivative, where the kernel is not acting on the concentration "
                       "variable directly.");
  params.addParam<Real>("gamma", 1.0, "Rate factor");
  return params;
}

RadialGreensSource::RadialGreensSource(const InputParameters & parameters)
  : Kernel(parameters),
    _convolution(getUserObject<RadialGreensConvolution>("convolution").getConvolution()),
    _is_coupled(isCoupled("c")),
    _c_var(_is_coupled ? coupled("c") : _var.number()),
    _c(_is_coupled ? coupledValue("c") : _u),
    _gamma(getParam<Real>("gamma")),
    _dt(_fe_problem.dt())
{
}

void
RadialGreensSource::precalculateResidual()
{
  _result = _convolution.find(_current_elem->id());
  if (_result == _convolution.end())
    mooseError("Current element not found in convolution result");
}

Real
RadialGreensSource::computeQpResidual()
{
  return _gamma * (_c[_qp] - _result->second[_qp] / (_JxW[_qp] * _coord[_qp])) * _test[_i][_qp];
}

Real
RadialGreensSource::computeQpJacobian()
{
  return _is_coupled ? 0.0 : computeQpCJacobian();
}

Real
RadialGreensSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  // c Off-Diagonal Jacobian
  if (_c_var == jvar)
    return computeQpCJacobian();

  return 0.0;
}

Real
RadialGreensSource::computeQpCJacobian()
{
  // Calculate the Jacobian for the c variable
  return _gamma * _phi[_j][_qp] * _test[_i][_qp];
}
