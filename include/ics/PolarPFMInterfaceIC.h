/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "InitialCondition.h"

/**
 * Upsilon IC, eq. (50) from Physical Review B 89, 184102 (2014)
 */
class PolarPFMInterfaceIC : public InitialCondition
{
public:
  static InputParameters validParams();

  PolarPFMInterfaceIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  /// coupled solid phase order parameter
  const VariableValue & _theta;

  /// interpolation function parameter
  const Real _a_beta;

  /// solid1 / melt gradient energy coefficient
  const Real _beta10;
  /// solid2 / melt gradient energy coefficient
  const Real _beta20;

  /// scale factor for interface width
  const Real _p;

  /// Minimum x coordinate
  const Real _xmin;

  /// Maximum x coordinate
  const Real _xmax;
};
