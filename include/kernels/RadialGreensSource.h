/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef RADIALGREENSSOURCE_H
#define RADIALGREENSSOURCE_H

#include "Kernel.h"
#include "RadialGreensConvolution.h"

class RadialGreensSource;

template <>
InputParameters validParams<RadialGreensSource>();

/**
 * Apply the convolution from a RadialGreensConvolution object to a non-linear variable
 */
class RadialGreensSource : public Kernel
{
public:
  RadialGreensSource(const InputParameters & parameters);

protected:
  void precalculateResidual() override;
  virtual Real computeQpResidual() override;

  // convolution result
  const RadialGreensConvolution::Result & _convolution;

  // rate factor
  const Real _gamma;

  // iterator pointing to the map entry for the current element
  RadialGreensConvolution::Result::const_iterator _result;
};

#endif // RADIALGREENSSOURCE_H
