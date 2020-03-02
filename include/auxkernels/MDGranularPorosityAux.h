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

class MDGranularPorosityAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MDGranularPorosityAux(const InputParameters & params);

  virtual Real computeValue();

protected:
  /// referene to the MDRunBase user object
  const MDRunBase & _md_uo;

  const bool _compute_packing;

  /// property value that is computed only on qp = 0
  Real _packing_fraction;
};
