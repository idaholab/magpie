/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMElementResultAccess.h"
#include "AuxKernel.h"

class MyTRIMElementResultAux : public MyTRIMElementResultAccess<AuxKernel>
{
public:
  static InputParameters validParams();

  MyTRIMElementResultAux(const InputParameters & params);

  virtual Real computeValue();
};
