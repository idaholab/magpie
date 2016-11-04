#include "AtomicDensityAux.h"
#include "MyTRIMRasterizer.h"

template<>
InputParameters validParams<AtomicDensityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Compute the atomic density at an element using data from a MyTRIMRasterizer");
  params.addRequiredParam<UserObjectName>("rasterizer", "MyTRIMRasterizer object to provide material data");
  return params;
}

AtomicDensityAux::AtomicDensityAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer"))
{
}

Real
AtomicDensityAux::computeValue()
{
  // prepare the material using data from the rasterizer
  const std::vector<Real> & material_data = _rasterizer.material(_current_elem);

  // calculate the lattice site occupation
  Real total_occupation = 0.0;
  for (auto t : material_data)
    total_occupation += t;

  // compute atomic density
  return total_occupation / _rasterizer.siteVolume(_current_elem);
}
