/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "TemporalCorrelation.h"

registerMooseObject("MagpieApp", TemporalCorrelation);

InputParameters
TemporalCorrelation::validParams()
{
  InputParameters params = ElementAverageValue::validParams();
  params.addClassDescription("Integrate over the L2 norm of a variable time derivative.");
  return params;
}

TemporalCorrelation::TemporalCorrelation(const InputParameters & parameters)
  : ElementAverageValue(parameters), _u_dot(MooseVariableInterface<Real>::dot())
{
}

Real
TemporalCorrelation::computeQpIntegral()
{
  return _u_dot[_qp] * _u_dot[_qp];
}

Real
TemporalCorrelation::getValue() const
{
  return std::sqrt(ElementAverageValue::getValue());
}
