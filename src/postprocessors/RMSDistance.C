/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "RMSDistance.h"

registerMooseObject("MooseApp", RMSDistance);

InputParameters
RMSDistance::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addRequiredParam<Point>("point", "Point to which the RMS distance is computed");
  params.addClassDescription(
      "Computes the RMS distance of a distribution with respect to a fixed point.");
  return params;
}

RMSDistance::RMSDistance(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters),
    _var(_fe_problem.getStandardVariable(_tid, getParam<std::vector<VariableName>>("variable")[0])),
    _point(getParam<Point>("point"))
{
}

Real
RMSDistance::computeQpIntegral()
{
  _normalization += _JxW[_qp] * _coord[_qp] * _u[_qp];
  return _mesh.minPeriodicVector(_var.number(), _point, _q_point[_qp]).norm_sq() * _u[_qp];
}

void
RMSDistance::initialize()
{
  ElementIntegralVariablePostprocessor::initialize();
  _normalization = 0;
}

Real
RMSDistance::getValue()
{
  Real integral = ElementIntegralVariablePostprocessor::getValue();

  gatherSum(_normalization);

  return std::sqrt(integral / _normalization);
}

void
RMSDistance::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  const RMSDistance & pps = static_cast<const RMSDistance &>(y);
  _normalization += pps._normalization;
}
