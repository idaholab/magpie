/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMPKAStatistics.h"
#include "MyTRIMRasterizer.h"
#include "MagpieParallel.h"
#include "mytrim/ion.h"

template<>
InputParameters validParams<MyTRIMPKAStatistics>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addClassDescription("Compile a table of the number of PKA with the same value of a chosen property (e.g. a mass histogram)");
  params.addRequiredParam<UserObjectName>("rasterizer", "Name of the MyTRIMRasterizer userobject to pull data from");
  MooseEnum value_type_options("MASS=0 ZAID");
  params.addParam<MooseEnum>("value_type", value_type_options, "The common property to bin the PKA set according to");
  return params;
}

MyTRIMPKAStatistics::MyTRIMPKAStatistics(const InputParameters & params) :
    GeneralVectorPostprocessor(params),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _value_type(getParam<MooseEnum>("value_type").getEnum<ValueType>()),
    _property(declareVector(getParam<MooseEnum>("value_type"))),
    _count(declareVector("n"))
{
}

void
MyTRIMPKAStatistics::initialize()
{
  _count_map.clear();
  _property.clear();
  _count.clear();
}

void
MyTRIMPKAStatistics::execute()
{
  const std::vector<MyTRIM_NS::IonBase> & pka_list = _rasterizer.getPKAList();

  unsigned int prop;
  for (auto & pka : pka_list)
  {
    const unsigned int Z = pka._Z;
    const unsigned int M = std::round(pka._m);;

    switch (_value_type)
    {
      case MASS:
        prop = M;
        break;

      case ZAID:
        prop = Z * 1000 + M;
        break;

      default:
        mooseError("Internal error");
    }

    auto i = _count_map.find(prop);
    if (i == _count_map.end())
      i = _count_map.insert(_count_map.begin(), std::make_pair(prop, 1));
    else
      i->second++;
  }
}

void
MyTRIMPKAStatistics::finalize()
{
  // for single processor runs we do not need to do anything here
  if (_app.n_processors() > 1)
  {
    // flatten the map
    std::vector<std::pair<unsigned int, Real>> count_flat(_count_map.begin(), _count_map.end());

    // allgather the flattend maps
    _communicator.allgather(count_flat);

    // reconstitute the map
    _count_map.clear();
    for (auto & count_pair: count_flat)
    {
      auto i = _count_map.find(count_pair.first);
      if (i == _count_map.end())
        i = _count_map.insert(_count_map.begin(), count_pair);
      else
        i->second += count_pair.second;
    }
  }

  // copy count map into count vector (ordering according to key is guaranteed)
  _property.reserve(_count_map.size());
  _count.reserve(_count_map.size());
  for (auto & count_pair: _count_map)
  {
    _property.push_back(count_pair.first);
    _count.push_back(count_pair.second);
  }
}
