/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMElementResultAux.h"
#include "MyTRIMElementRun.h"

registerMooseObject("MagpieApp", MyTRIMElementResultAux);

template<>
InputParameters validParams<MyTRIMElementResultAux>()
{
  InputParameters params = MyTRIMElementResultAccess<AuxKernel>::validParams();
  return params;
}

MyTRIMElementResultAux::MyTRIMElementResultAux(const InputParameters & parameters) :
    MyTRIMElementResultAccess<AuxKernel>(parameters)
{
  /**
   * having this AuxKernel also depend on the rasterizer bumps the rasterizer into
   * the preaux group and ensures it is executed _before_ the MyTRIMRun object.
   */
  getUserObjectByName<MyTRIMRasterizer>(_mytrim.getRasterizerName());
}

Real
MyTRIMElementResultAux::computeValue()
{
  return getDefectRate();
}
