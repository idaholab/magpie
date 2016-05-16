#include "MyTRIMElementSource.h"
#include "MyTRIMElementRun.h"

template<>
InputParameters validParams<MyTRIMElementSource>()
{
  InputParameters params = MyTRIMResultAccess<Kernel>::validParams();
  params.addParam<Real>("prefactor", 1.0, "Prefactor to scale the applied production rate");
  return params;
}

MyTRIMElementSource::MyTRIMElementSource(const InputParameters & parameters) :
    MyTRIMResultAccess<Kernel>(parameters),
    _prefactor(getParam<Real>("prefactor"))
{
}

Real
MyTRIMElementSource::computeQpResidual()
{
  return -_prefactor * getDefectRate() * _test[_i][_qp];
}
