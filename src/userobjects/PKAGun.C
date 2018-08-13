/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAGun.h"

registerMooseObject("MagpieApp", PKAGun);

template <>
InputParameters
validParams<PKAGun>()
{
  InputParameters params = validParams<PKAFixedPointGenerator>();
  params.addClassDescription(
      "This PKAGenerator starts particle from a fixed point in a fixed direction.");
  params.addRequiredParam<Point>("direction", "The fixed direction the PKAs move along");
  return params;
}

PKAGun::PKAGun(const InputParameters & parameters)
  : PKAFixedPointGenerator(parameters), _direction(getParam<Point>("direction"))
{
}

void
PKAGun::setDirection(MyTRIM_NS::IonBase & ion) const
{
  ion._dir = _direction;
}
