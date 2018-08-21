/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAList.h"
#include "MyTRIMRasterizer.h"

registerMooseObject("MagpieApp", PKAList);

template <>
InputParameters
validParams<PKAList>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addClassDescription("Dumps the entire PKA list");
  params.addRequiredParam<UserObjectName>(
      "rasterizer", "Name of the MyTRIMRasterizer userobject that provides the PKA list.");
  return params;
}

PKAList::PKAList(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _pka_list(_rasterizer.getPKAList()),
    _x(declareVector("x")),
    _y(declareVector("y")),
    _z(declareVector("z")),
    _seed(declareVector("seed")),
    _m(declareVector("m")),
    _Z(declareVector("Z"))
{
}

void
PKAList::initialize()
{
  _x.clear();
  _y.clear();
  _z.clear();
  _seed.clear();
  _m.clear();
  _Z.clear();
}

void
PKAList::execute()
{
  for (auto & pka : _pka_list)
  {
    _x.push_back(pka._pos(0));
    _y.push_back(pka._pos(1));
    _z.push_back(pka._pos(2));
    _seed.push_back(pka._seed);
    _m.push_back(pka._m);
    _Z.push_back(pka._Z);
  }
}

void
PKAList::finalize()
{
  // broadcast data to processor 0
  _communicator.allgather(_x, false);
  _communicator.allgather(_y, false);
  _communicator.allgather(_z, false);
  _communicator.allgather(_seed, false);
  _communicator.allgather(_m, false);
  _communicator.allgather(_Z, false);
}
