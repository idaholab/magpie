/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include <MooseFFTVariable.h>

InputParameters
MooseFFTVariable::validParams()
{
  auto params = MooseVariableFEBase::validParams();
  return params;
}

MooseFFTVariable::MooseFFTVariable(const InputParameters & parameters) : MooseVariable(parameters)
{
}
