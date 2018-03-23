/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMDiracSource.h"
#include "MyTRIMDiracRun.h"
#include "MooseMesh.h"

registerMooseObject("MagpieApp", MyTRIMDiracSource);

template<>
InputParameters validParams<MyTRIMDiracSource>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMDiracRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");

  // this needs multiplicity enabled to give meaningful results
  params.set<bool>("drop_duplicate_points") = false;

  return params;
}

MyTRIMDiracSource::MyTRIMDiracSource(const InputParameters & parameters) :
    DiracKernel(parameters),
    _mytrim(getUserObject<MyTRIMDiracRun>("runner")),
    _ivar(getParam<unsigned int>("ivar")),
    _defect(getParam<MooseEnum>("defect"))
{
  if (getParam<bool>("drop_duplicate_points") == true)
    mooseWarning("Explicitly setting drop_duplicate_points to true will cause overlapping defects to be miscounted.");
}

void
MyTRIMDiracSource::addPoints()
{
  for (auto && defect : _mytrim.result())
    if (defect._type == _defect && defect._var == _ivar)
    {
      Elem * elem = _mesh.elemPtr(defect._elem_id);
      if (elem)
        addPoint(elem, defect._location);
    }
}

Real
MyTRIMDiracSource::computeQpResidual()
{
  return -_test[_i][_qp];
}
