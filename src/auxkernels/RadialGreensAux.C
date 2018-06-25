/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "RadialGreensAux.h"
#include "RadialGreensConvolution.h"

registerMooseObject("MagpieApp", RadialGreensAux);

template <>
InputParameters
validParams<RadialGreensAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Visualize data generated in a RadialGreensConvolution user object");
  params.addRequiredParam<UserObjectName>("convolution", "RadialGreensConvolution user object");
  return params;
}

RadialGreensAux::RadialGreensAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _convolution(getUserObject<RadialGreensConvolution>("convolution").getConvolution())
{
  if (isNodal())
    paramError("variable", "RadialGreensAux must be applied to an elemental AuxVariable");
}

Real
RadialGreensAux::computeValue()
{
  auto it = _convolution.find(_current_elem->id());
  if (it != _convolution.end())
    return it->second[_qp] / (_JxW[_qp] * _coord[_qp]);

  return 0.0;
}
