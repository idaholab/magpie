#include "MyTRIMElementHeatSource.h"
#include "MyTRIMElementRun.h"

template<>
InputParameters validParams<MyTRIMElementHeatSource>()
{
  InputParameters params = MyTRIMElementEnergyAccess<Kernel>::validParams();
  return params;
}

MyTRIMElementHeatSource::MyTRIMElementHeatSource(const InputParameters & parameters) :
    MyTRIMElementEnergyAccess<Kernel>(parameters)
{
}

Real
MyTRIMElementHeatSource::computeQpResidual()
{
  return -getEnergyDensity() * _test[_i][_qp];
}
