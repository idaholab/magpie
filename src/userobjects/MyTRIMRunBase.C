/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMRunBase.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<MyTRIMRunBase>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<UserObjectName>("rasterizer", "MyTRIMRasterizer object to provide material data");

  // we run this object once a timestep
  params.set<MultiMooseEnum>("execute_on") = "timestep_begin";
  params.suppressParameter<MultiMooseEnum>("execute_on");

  return params;
}

MyTRIMRunBase::MyTRIMRunBase(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _nvars(_rasterizer.nVars()),
    _pka_list(_rasterizer.getPKAList()),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension())
{
  if (_mesh.isParallelMesh())
    mooseError("MyTRIM runs currently require serial meshes.");

  if (_dim < 2 || _dim > 3)
    mooseError("TRIM simulation works in 2D or 3D only.");

  if (_dim == 2)
  {
    // make sure all nodes lie in the xy-plane
    const auto nd_end = _mesh.getMesh().nodes_end();
    for (auto nd = _mesh.getMesh().nodes_begin(); nd != nd_end; ++nd)
      if ((**nd)(2) != 0.0)
        mooseError("Two dimensional meshes must lie in the z=0 plane.");
  }
}
