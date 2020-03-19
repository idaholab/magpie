/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"

/**
 * Get a selected component of the gradient of a coupled variable
 */
class GradientComponentAux : public AuxKernel
{
public:
  static InputParameters validParams();

  GradientComponentAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Gradient of the coupled variable
  const VariableGradient & _grad_v;

  /// Component of the gradient vector to match
  const unsigned int _component;
};