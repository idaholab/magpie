/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMMesh.h"
#include "PointLocatorRegularOrthogonal.h"

registerMooseObject("MagpieApp", MyTRIMMesh);

InputParameters
MyTRIMMesh::validParams()
{
  InputParameters params = GeneratedMesh::validParams();
  params.addClassDescription(
      "Regular orthogonal generated mesh with restrictions to simple elements.");
  params.suppressParameter<Real>("bias_x");
  params.suppressParameter<Real>("bias_y");
  params.suppressParameter<Real>("bias_z");
  return params;
}

MyTRIMMesh::MyTRIMMesh(const InputParameters & parameters)
  : GeneratedMesh(parameters),
    _cell_count({_nx, _ny, _nz}),
    _min_corner(_xmin, _ymin, _zmin),
    _max_corner(_xmax, _ymax, _zmax)
{
  _cell_count.resize(_dim);

  if (isParamValid("elem_type"))
  {
    MooseEnum elem_type_enum = getParam<MooseEnum>("elem_type");

    // validate element type (only orthogonal stuff)
    switch (_dim)
    {
      case 1:
        if (elem_type_enum != "EDGE2")
          mooseError("1D simulations need to use EDGE2 elements");
        break;
      case 2:
        if (elem_type_enum != "QUAD4")
          mooseError("2D simulations need to use QUAD4 elements");
        break;
      case 3:
        if (elem_type_enum != "HEX8")
          mooseError("3D simulations need to use HEX8 elements");
    }
  }

  // actually let's not allow 1D for now
  if (_dim == 1)
    mooseError("Only 2D and 3D simulations are currently supported.");

  // make sure no biasing snuck in
  if (_bias_x != 1.0 || _bias_y != 1.0 || _bias_z != 1.0)
    mooseError("Biased meshes are not supported.");
}

MyTRIMMesh::MyTRIMMesh(const MyTRIMMesh & other_mesh)
  : GeneratedMesh(other_mesh), _cell_count(other_mesh._cell_count)
{
}

std::unique_ptr<MooseMesh>
MyTRIMMesh::safeClone() const
{
  return std::make_unique<MyTRIMMesh>(*this);
}

std::unique_ptr<PointLocatorBase>
MyTRIMMesh::getPointLocator() const
{
  if (_point_locator.get() == nullptr)
  {
    _point_locator.reset(new PointLocatorRegularOrthogonal(getMesh()));
    _point_locator->init(_cell_count, _min_corner, _max_corner);
  }
  return std::unique_ptr<PointLocatorBase>(
      new PointLocatorRegularOrthogonal(getMesh(), _point_locator.get()));
}

unsigned int
MyTRIMMesh::getCellCountInDimension(unsigned int component)
{
  mooseAssert(component < dimension(),
              "Querying invalid component in MyTRIMMesh::getCellCountInDimension");
  return _cell_count[component];
}
