/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "Action.h"

/**
 * Sets up the polar phase field model by Momeni et al.
 */
class PolarPhaseFieldAction : public Action
{
public:
  static InputParameters validParams();

  PolarPhaseFieldAction(const InputParameters & parameters);

  virtual void act() override;
};
