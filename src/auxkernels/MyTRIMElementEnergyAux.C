#include "MyTRIMElementEnergyAux.h"
#include "MyTRIMElementRun.h"

template<>
InputParameters validParams<MyTRIMElementEnergyAux>()
{
  InputParameters params = MyTRIMElementEnergyAccess<AuxKernel>::validParams();
  return params;
}

MyTRIMElementEnergyAux::MyTRIMElementEnergyAux(const InputParameters & parameters) :
    MyTRIMElementEnergyAccess<AuxKernel>(parameters)
{
  /**
   * having this AuxKernel also depend on the rasterizer bumps the rasterizer into
   * the preaux group and ensures it is executed _before_ the MyTRIMRun object.
   */
  getUserObjectByName<MyTRIMRasterizer>(_mytrim.getRasterizerName());
}

Real
MyTRIMElementEnergyAux::computeValue()
{
  return getEnergyDensity();
}
