/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "AtomicDensityAux.h"
#include "MyTRIMRasterizer.h"

#include "libmesh/utility.h"

registerMooseObject("MagpieApp", AtomicDensityAux);

template <>
InputParameters
validParams<AtomicDensityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Compute the atomic density as Atoms/volume at an element using data "
                             "from a MyTRIMRasterizer. Volume is in the Mesh units set in the "
                             "rasterizer.");
  params.addRequiredParam<UserObjectName>("rasterizer",
                                          "MyTRIMRasterizer object to provide material data");
  return params;
}

AtomicDensityAux::AtomicDensityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _volume_scale(Utility::pow<3>(_rasterizer.getTrimParameters().length_scale / 1000.0))
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

  // compute atomic density (siteVolume is in nm^3)
  return _volume_scale * total_occupation / _rasterizer.siteVolume(_current_elem);
}
