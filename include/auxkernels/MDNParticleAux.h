/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"

class MDRunBase;

class MDNParticleAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MDNParticleAux(const InputParameters & params);

  virtual Real computeValue();

protected:
  const MDRunBase & _md_uo;

  std::vector<unsigned int> _particles;
};
