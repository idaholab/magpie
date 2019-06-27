/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"

// forward declarations
class SPPARKSUserObject;
class SPPARKSAux;

template <>
InputParameters validParams<SPPARKSAux>();

class SPPARKSAux : public AuxKernel
{
public:
  SPPARKSAux(const InputParameters & params);
  virtual ~SPPARKSAux() {}

  virtual Real computeValue();

protected:
  const SPPARKSUserObject & _spparks;
  const unsigned int _var;
  const unsigned int _array;
};

