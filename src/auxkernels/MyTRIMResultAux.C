#include "MyTRIMResultAux.h"
#include "MyTRIMRun.h"

template<>
InputParameters validParams<MyTRIMResultAux>()
{
  InputParameters params = MyTRIMResultAccess<AuxKernel>::validParams();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMRun userobject to pull data from.");
  return params;
}

MyTRIMResultAux::MyTRIMResultAux(const InputParameters & parameters) :
    MyTRIMResultAccess<AuxKernel>(parameters)
{
}

Real
MyTRIMResultAux::computeValue()
{
  return getDefectRate();
}
