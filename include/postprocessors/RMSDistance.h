/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ElementIntegralVariablePostprocessor.h"

/**
 * This postprocessor computes The RMS distance of a variable to a given point
 */
class RMSDistance : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  RMSDistance(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

  using Postprocessor::getValue;
  virtual Real getValue() const override;

protected:
  virtual Real computeQpIntegral() override;

  MooseVariable & _var;
  const Point _point;
  Real _normalization;
};
