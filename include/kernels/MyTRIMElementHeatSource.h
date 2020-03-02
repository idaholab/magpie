/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMElementEnergyAccess.h"
#include "Kernel.h"

class MyTRIMElementHeatSource : public MyTRIMElementEnergyAccess<Kernel>
{
public:
  static InputParameters validParams();

  MyTRIMElementHeatSource(const InputParameters & params);

protected:
  virtual Real computeQpResidual();

  // current timestep size
  const Real & _dt;
};
