#include "MyTRIMSource.h"
#include "MyTRIMRun.h"

template<>
InputParameters validParams<MyTRIMSource>()
{
  InputParameters params = MyTRIMResultAccess<Kernel>::validParams();
  params.addParam<Real>("prefactor", 1.0, "Prefactor to scale the applied production rate");
  return params;
}

MyTRIMSource::MyTRIMSource(const InputParameters & parameters) :
    MyTRIMResultAccess<Kernel>(parameters),
    _prefactor(getParam<Real>("prefactor"))
{
}

Real
MyTRIMSource::computeQpResidual()
{
  return -_prefactor * getDefectRate() * _test[_i][_qp];
}
