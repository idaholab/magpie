/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"
#include "RadialGreensConvolution.h"

/**
 * Visualize data generated in a RadialGreensConvolution user object
 */
class RadialGreensAux : public AuxKernel
{
public:
  static InputParameters validParams();

  RadialGreensAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const RadialGreensConvolution::Result & _convolution;
};
