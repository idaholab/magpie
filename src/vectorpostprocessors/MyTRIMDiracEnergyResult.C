/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMDiracEnergyResult.h"
#include "MyTRIMDiracRun.h"

registerMooseObject("MagpieApp", MyTRIMDiracEnergyResult);

template <>
InputParameters
validParams<MyTRIMDiracEnergyResult>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addRequiredParam<UserObjectName>(
      "runner", "Name of the MyTRIMDiracRun userobject to pull data from.");
  return params;
}

MyTRIMDiracEnergyResult::MyTRIMDiracEnergyResult(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _mytrim(getUserObject<MyTRIMDiracRun>("runner")),
    _x(declareVector("x")),
    _y(declareVector("y")),
    _z(declareVector("z")),
    _elem_id(declareVector("elem_id")),
    _energy_deposition(declareVector("energy_deposition"))
{
}

void
MyTRIMDiracEnergyResult::initialize()
{
  _x.clear();
  _y.clear();
  _z.clear();
  _elem_id.clear();
  _energy_deposition.clear();
}

void
MyTRIMDiracEnergyResult::execute()
{
  for (auto && entry : _mytrim.result())
    if (static_cast<unsigned int>(entry._type) == 5)
    {
      _x.push_back(entry._location(0));
      _y.push_back(entry._location(1));
      _z.push_back(entry._location(2));
      _elem_id.push_back(entry._elem_id);
      _energy_deposition.push_back(entry._weight);
    }
}

void
MyTRIMDiracEnergyResult::finalize()
{
  // the MyTRIMDiracRunner already does the necessary communication
}
