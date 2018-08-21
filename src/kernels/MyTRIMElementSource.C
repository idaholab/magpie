/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMElementSource.h"
#include "MyTRIMElementRun.h"

registerMooseObject("MagpieApp", MyTRIMElementSource);

template <>
InputParameters
validParams<MyTRIMElementSource>()
{
  InputParameters params = MyTRIMElementResultAccess<Kernel>::validParams();
  params.addParam<Real>("prefactor", 1.0, "Prefactor to scale the applied production rate");
  return params;
}

MyTRIMElementSource::MyTRIMElementSource(const InputParameters & parameters)
  : MyTRIMElementResultAccess<Kernel>(parameters),
    _prefactor(getParam<Real>("prefactor")),
    _trim_parameters(_rasterizer.getTrimParameters())
{
}

Real
MyTRIMElementSource::computeQpResidual()
{
  return -_prefactor * getDefectRate() / _trim_parameters.last_executed_dt * _test[_i][_qp];
}
