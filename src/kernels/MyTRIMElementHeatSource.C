/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMElementHeatSource.h"
#include "MyTRIMElementRun.h"

registerMooseObject("MagpieApp", MyTRIMElementHeatSource);

template <>
InputParameters
validParams<MyTRIMElementHeatSource>()
{
  InputParameters params = MyTRIMElementEnergyAccess<Kernel>::validParams();
  return params;
}

MyTRIMElementHeatSource::MyTRIMElementHeatSource(const InputParameters & parameters)
  : MyTRIMElementEnergyAccess<Kernel>(parameters), _dt(_fe_problem.dt())
{
}

Real
MyTRIMElementHeatSource::computeQpResidual()
{
  return -getEnergyDensity() / _dt * _test[_i][_qp];
}
