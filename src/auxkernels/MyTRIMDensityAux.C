/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMDensityAux.h"
#include "MooseMyTRIMMaterial.h"

registerMooseObject("MagpieApp", MyTRIMDensityAux);

template <>
InputParameters
validParams<MyTRIMDensityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Returns the material density in g/cm^3");
  params.addRequiredParam<UserObjectName>("rasterizer",
                                          "MyTRIMRasterizer object to provide material data");
  return params;
}

MyTRIMDensityAux::MyTRIMDensityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _trim_parameters(_rasterizer.getTrimParameters()),
    _nvars(_trim_parameters.nVars())
{
  if (isNodal())
    mooseError("MyTRIMDensityAux needs to be applied to an elemental AuxVariable.");
}

Real
MyTRIMDensityAux::computeValue()
{
  if (_qp == 0)
  {
    // prepare the material using data from the rasterizer
    const std::vector<Real> & material_data = _rasterizer.material(_current_elem);
    MooseMyTRIMMaterial material(&_simconf);

    // set elements
    MyTRIM_NS::Element element;
    for (unsigned int i = 0; i < _nvars; ++i)
    {
      element = _trim_parameters.element_prototypes[i];
      element._t = material_data[i];
      material._element.push_back(element);
    }

    // calculate the density (must happen first!)
    material.calculateDensity(_rasterizer.siteVolume(_current_elem));

    _value_cache = material._rho;
  }

  return _value_cache;
}
