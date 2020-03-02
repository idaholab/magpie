/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMDiracResult.h"
#include "MyTRIMDiracRun.h"

registerMooseObject("MagpieApp", MyTRIMDiracResult);

InputParameters
MyTRIMDiracResult::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>(
      "runner", "Name of the MyTRIMDiracRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");
  return params;
}

MyTRIMDiracResult::MyTRIMDiracResult(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _mytrim(getUserObject<MyTRIMDiracRun>("runner")),
    _ivar(getParam<unsigned int>("ivar")),
    _defect(getParam<MooseEnum>("defect")),
    _x(declareVector("x")),
    _y(declareVector("y")),
    _z(declareVector("z")),
    _elem_id(declareVector("elem_id"))
{
}

void
MyTRIMDiracResult::initialize()
{
  _x.clear();
  _y.clear();
  _z.clear();
  _elem_id.clear();
}

void
MyTRIMDiracResult::execute()
{
  for (auto && defect : _mytrim.result())
    if (defect._type >= 0 && static_cast<unsigned int>(defect._type) == _defect &&
        defect._var == _ivar)
    {
      _x.push_back(defect._location(0));
      _y.push_back(defect._location(1));
      _z.push_back(defect._location(2));
      _elem_id.push_back(defect._elem_id);
    }
}

void
MyTRIMDiracResult::finalize()
{
  // the MyTRIMDiracRunner already does the necessary communication
}
