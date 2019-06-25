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
class MDRunBase;
class MDNParticleAux;

template <>
InputParameters validParams<MDNParticleAux>();

class MDNParticleAux : public AuxKernel
{
public:
  MDNParticleAux(const InputParameters & params);
  virtual ~MDNParticleAux() {}

  virtual Real computeValue();

protected:
  const MDRunBase & _md_uo;
};

