#include "MyTRIMElementResultAux.h"
#include "MyTRIMElementRun.h"

template<>
InputParameters validParams<MyTRIMElementResultAux>()
{
  InputParameters params = MyTRIMResultAccess<AuxKernel>::validParams();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMElementRun userobject to pull data from.");
  return params;
}

MyTRIMElementResultAux::MyTRIMElementResultAux(const InputParameters & parameters) :
    MyTRIMResultAccess<AuxKernel>(parameters)
{
}

Real
MyTRIMElementResultAux::computeValue()
{
  return getDefectRate();
}
