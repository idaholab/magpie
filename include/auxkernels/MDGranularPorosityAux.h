/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MDGRANULARPOROSITYAUX_H
#define MDGRANULARPOROSITYAUX_H

#include "AuxKernel.h"

// forward declarations
class MDRunBase;
class MDGranularPorosityAux;

template <>
InputParameters validParams<MDGranularPorosityAux>();

class MDGranularPorosityAux : public AuxKernel
{
public:
  MDGranularPorosityAux(const InputParameters & params);
  virtual ~MDGranularPorosityAux() {}

  virtual Real computeValue();

protected:
  /// referene to the MDRunBase user object
  const MDRunBase & _md_uo;

  const bool _compute_packing;

  /// property value that is computed only on qp = 0
  Real _packing_fraction;
};

#endif // MDGRANULARPOROSITYAUXx_H
