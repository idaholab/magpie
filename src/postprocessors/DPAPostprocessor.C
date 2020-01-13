/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef GSL_ENABLED

#include "DPAPostprocessor.h"
#include "DPAUserObjectBase.h"

registerMooseObject("MagpieApp", DPAPostprocessor);

template <>
InputParameters
validParams<DPAPostprocessor>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<UserObjectName>("dpa_object", "The neutronics damage object.");
  params.addClassDescription("Retrieves the value of the dpa from a DPAUserObjectBase.");
  return params;
}

DPAPostprocessor::DPAPostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params), _damage_object(getUserObject<DPAUserObjectBase>("dpa_object"))
{
}

Real
DPAPostprocessor::getValue()
{
  return _damage_object.getDPA();
}

#endif
