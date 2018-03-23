/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "SPPARKSAux.h"
#include "SPPARKSUserObject.h"

registerMooseObject("MagpieApp", SPPARKSAux);

template<>
InputParameters validParams<SPPARKSAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("user_object", "Name of SPPARKSUserObject");
  params.addParam<unsigned int>("var", "Index into SPPARKS array");
  MooseEnum arrayType("IARRAY DARRAY");
  params.addParam<MooseEnum>("array", arrayType, "SPPARKS array to read from");
  return params;
}

SPPARKSAux::SPPARKSAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _spparks(getUserObject<SPPARKSUserObject>("user_object")),
    _var(getParam<unsigned int>("var")),
    _array(getParam<MooseEnum>("array"))
{
}

Real
SPPARKSAux::computeValue()
{
  switch (_array)
  {
    case 0: // Integer Array
      return  _spparks.getIntValue(_current_node->id(), _var);

    case 1: // Double Array
      return  _spparks.getDoubleValue(_current_node->id(), _var);

    default:
      mooseError("Internal error.");
  }
}
