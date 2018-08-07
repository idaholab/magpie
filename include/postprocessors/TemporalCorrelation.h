/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef TEMPORALCORRELATION_H
#define TEMPORALCORRELATION_H

#include "ElementAverageValue.h"

class TemporalCorrelation;

template <>
InputParameters validParams<TemporalCorrelation>();

/**
 * Compute average of the L2 norm of a variable time derivative as a measure of the
 * correlation between timesteps. Small values indicate slow, gradual change, large
 * values indicate strong fluctuations.
 */
class TemporalCorrelation : public ElementAverageValue
{
public:
  TemporalCorrelation(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;
  virtual Real getValue() override;
};

#endif // TEMPORALCORRELATION_H
