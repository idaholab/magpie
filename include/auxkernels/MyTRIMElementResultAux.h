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

// forward declarations
class MyTRIMElementRun;
class MyTRIMElementResultAux;

template <>
InputParameters validParams<MyTRIMElementResultAux>();

class MyTRIMElementResultAux : public MyTRIMElementResultAccess<AuxKernel>
{
public:
  MyTRIMElementResultAux(const InputParameters & params);
  virtual ~MyTRIMElementResultAux() {}

  virtual Real computeValue();
};

