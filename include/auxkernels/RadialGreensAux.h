/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef RADIALGREENSAUX_H
#define RADIALGREENSAUX_H

#include "AuxKernel.h"
#include "RadialGreensConvolution.h"

class RadialGreensAux;

template <>
InputParameters validParams<RadialGreensAux>();

/**
 * Visualize data generated in a RadialGreensConvolution user object
 */
class RadialGreensAux : public AuxKernel
{
public:
  RadialGreensAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const RadialGreensConvolution::Result & _convolution;
};

#endif // RADIALGREENSAUX_H
