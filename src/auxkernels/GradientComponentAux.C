/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "GradientComponentAux.h"

registerMooseObject("MagpieApp", GradientComponentAux);

InputParameters
GradientComponentAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Return specified component of the gradient of a coupled variable.");
  params.addRequiredCoupledVar("v", "Coupled variable to extract gradient component of");

  params.addRequiredParam<unsigned int>("component",
                                        "Component of the gradient of the coupled variable v");
  return params;
}

GradientComponentAux::GradientComponentAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_v(coupledGradient("v")),
    _component(getParam<unsigned int>("component"))
{
  if (_component >= LIBMESH_DIM)
    paramError("component", "Component too large for LIBMESH_DIM");
}

Real
GradientComponentAux::computeValue()
{
  return _grad_v[_qp](_component);
}
