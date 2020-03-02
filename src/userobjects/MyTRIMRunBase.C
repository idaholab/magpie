/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMRunBase.h"
#include "MooseMesh.h"

InputParameters
MyTRIMRunBase::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addRequiredParam<UserObjectName>("rasterizer",
                                          "MyTRIMRasterizer object to provide material data");

  // we run this object once a timestep
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  params.suppressParameter<ExecFlagEnum>("execute_on");

  return params;
}

MyTRIMRunBase::MyTRIMRunBase(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _trim_parameters(_rasterizer.getTrimParameters()),
    _nvars(_trim_parameters.nVars()),
    _pka_list(_rasterizer.getPKAList()),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension())
{
  if (_mesh.isDistributedMesh())
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
