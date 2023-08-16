/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ElementAverageValue.h"

/**
 * Compute average of the L2 norm of a variable time derivative as a measure of the
 * correlation between timesteps. Small values indicate slow, gradual change, large
 * values indicate strong fluctuations.
 */
class TemporalCorrelation : public ElementAverageValue
{
public:
  static InputParameters validParams();

  TemporalCorrelation(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  using Postprocessor::getValue;
  virtual Real getValue() const override;

  const VariableValue & _u_dot;
};
