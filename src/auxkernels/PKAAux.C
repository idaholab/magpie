#include "PKAAux.h"
#include "MyTRIMRasterizer.h"

template<>
InputParameters validParams<PKAAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("rasterizer", "MyTRIMRasterizer object to provide material data");
  return params;
}

PKAAux::PKAAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer"))
{
}

Real
PKAAux::computeValue()
{
  // we have to divide by the volume to get the number and not number * volume
  return _rasterizer.getNPKA(_current_elem->id()) / _current_elem->volume();
}
