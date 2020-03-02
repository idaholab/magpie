/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPhaseFieldAction.h"

registerMooseAction("MagpieApp", PolarPhaseFieldAction, "add_kernel");

InputParameters
PolarPhaseFieldAction::validParams()
{
  InputParameters params = Action::validParams();
  return params;
}

PolarPhaseFieldAction::PolarPhaseFieldAction(const InputParameters & parameters)
  : Action(parameters)
{
}

void
PolarPhaseFieldAction::act()
{
}
