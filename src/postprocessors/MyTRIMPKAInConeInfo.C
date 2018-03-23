/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMPKAInConeInfo.h"

registerMooseObject("MagpieApp", MyTRIMPKAInConeInfo);

template<>
InputParameters validParams<MyTRIMPKAInConeInfo>()
{
  InputParameters params = validParams<MyTRIMPKAInfo>();
  params.addClassDescription("Aggregate a global property of the primary knock-on atom for PKAs within a cone of "
                             "opening angle opening_angle along cone_axis");
  params.addRequiredParam<RealVectorValue>("cone_axis", "Axis of the cone");
  params.addRequiredParam<Real>("opening_angle", "Opening angle of the cone = twice the angle measured from cone_axis to surface");
  return params;
}

MyTRIMPKAInConeInfo::MyTRIMPKAInConeInfo(const InputParameters & params) :
    MyTRIMPKAInfo(params),
    _direction(getParam<RealVectorValue>("cone_axis")),
    _min_cosine(std::cos(0.5 * getParam<Real>("opening_angle")))
{
}

bool
MyTRIMPKAInConeInfo::skipPKA(const MyTRIM_NS::IonBase & ion) const
{
  Real magnitude = _direction.norm();
  magnitude *= ion._dir.norm();
  // we skip counting if the cosine is smaller than min cosine because that
  // translates to an angle larger than opening_angle / 2
  return _direction * ion._dir < magnitude * _min_cosine;
}
