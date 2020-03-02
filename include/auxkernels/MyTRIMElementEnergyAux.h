/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMElementEnergyAccess.h"
#include "AuxKernel.h"

class MyTRIMElementEnergyAux : public MyTRIMElementEnergyAccess<AuxKernel>
{
public:
  static InputParameters validParams();

  MyTRIMElementEnergyAux(const InputParameters & params);

  virtual Real computeValue();
};
