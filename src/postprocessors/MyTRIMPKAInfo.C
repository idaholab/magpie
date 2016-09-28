#include "MyTRIMPKAInfo.h"
#include "MyTRIMRasterizer.h"
#include "MagpieParallel.h"
#include "mytrim/ion.h"

template<>
InputParameters validParams<MyTRIMPKAInfo>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addClassDescription("Aggregate a global property of the primary knock-on atom (PKA) list (e.g. total energy or number of PKA)");
  params.addRequiredParam<UserObjectName>("rasterizer", "Name of the MyTRIMRasterizer userobject to pull data from");
  MooseEnum value_type_options("TOTAL_MASS=0 TOTAL_ENERGY TOTAL_CHARGE TOTAL_NUMBER");
  params.addParam<MooseEnum>("value_type", value_type_options, "The property of the PKA set which is aggregated by this postprocessor");
  return params;
}

MyTRIMPKAInfo::MyTRIMPKAInfo(const InputParameters & params) :
    GeneralPostprocessor(params),
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
  switch (_value_type)
  {
    case TOTAL_MASS:
      for (auto & pka : pka_list)
        _value += pka._m;
      break;

    case TOTAL_ENERGY:
      for (auto & pka : pka_list)
        _value += pka._E;
      break;

    case TOTAL_CHARGE:
      for (auto & pka : pka_list)
        _value += pka._Z;
      break;

    case TOTAL_NUMBER:
      _value += pka_list.size();
      break;

    default:
      mooseError("Internal error");
  }
}

void
MyTRIMPKAInfo::finalize()
{
  gatherSum(_value);
}

Real
MyTRIMPKAInfo::getValue()
{
  return _value;
}
