/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMPKAInfo.h"
#include "MyTRIMRasterizer.h"

registerMooseObject("MagpieApp", MyTRIMPKAInfo);

InputParameters
MyTRIMPKAInfo::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription("Aggregate a global property of the primary knock-on atom (PKA) list "
                             "(e.g. total energy or number of PKA)");
  params.addRequiredParam<UserObjectName>(
      "rasterizer", "Name of the MyTRIMRasterizer userobject to pull data from");
  MooseEnum value_type_options("TOTAL_MASS=0 TOTAL_ENERGY TOTAL_CHARGE TOTAL_NUMBER");
  params.addParam<MooseEnum>(
      "value_type",
      value_type_options,
      "The property of the PKA set which is aggregated by this postprocessor");
  return params;
}

MyTRIMPKAInfo::MyTRIMPKAInfo(const InputParameters & params)
  : GeneralPostprocessor(params),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _value_type(getParam<MooseEnum>("value_type").getEnum<ValueType>())
{
}

void
MyTRIMPKAInfo::initialize()
{
  _value = 0.0;
}

void
MyTRIMPKAInfo::execute()
{
  const std::vector<MyTRIM_NS::IonBase> & pka_list = _rasterizer.getPKAList();
  for (auto & pka : pka_list)
  {
    if (skipPKA(pka))
      continue;

    switch (_value_type)
    {
      case TOTAL_MASS:
        _value += pka._m;
        break;

      case TOTAL_ENERGY:
        _value += pka._E;
        break;

      case TOTAL_CHARGE:
        _value += pka._Z;
        break;

      case TOTAL_NUMBER:
        _value += 1.0;
        break;

      default:
        mooseError("Internal error");
    }
  }
}

void
MyTRIMPKAInfo::finalize()
{
  gatherSum(_value);
}

Real
MyTRIMPKAInfo::getValue() const
{
  return _value;
}
