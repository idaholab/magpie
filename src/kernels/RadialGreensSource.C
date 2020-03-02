/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "RadialGreensSource.h"

registerMooseObject("MagpieApp", RadialGreensSource);

InputParameters
RadialGreensSource::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Apply the convolution from a RadialGreensConvolution object to a non-linear variable");
  params.addRequiredParam<UserObjectName>("convolution", "RadialGreensConvolution user object");
  params.addParam<Real>("gamma", 1.0, "Rate factor");
  return params;
}

RadialGreensSource::RadialGreensSource(const InputParameters & parameters)
  : Kernel(parameters),
    _convolution(getUserObject<RadialGreensConvolution>("convolution").getConvolution()),
    _gamma(getParam<Real>("gamma"))
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
  return -_gamma * _result->second[_qp] / (_JxW[_qp] * _coord[_qp]) * _test[_i][_qp];
}
