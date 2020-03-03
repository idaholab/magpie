/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef GSL_ENABLED

#pragma once

#include "GeneralPostprocessor.h"

// forward declarations
class DPAUserObjectBase;

class DPAPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  DPAPostprocessor(const InputParameters & parameters);
  virtual void execute() override {}
  virtual void initialize() override {}
  virtual Real getValue() override;

protected:
  const DPAUserObjectBase & _damage_object;
};

#endif
