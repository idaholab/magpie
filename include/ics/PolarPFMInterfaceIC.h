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
  const VariableValue & _theta;

  const Real _a_beta;
  const Real _beta10;
  const Real _beta20;
};
