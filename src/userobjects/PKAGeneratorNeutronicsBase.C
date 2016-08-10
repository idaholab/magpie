#include "PKAGeneratorNeutronicsBase.h"
#include "MagpieUtils.h"
#include "MultiIndex.h"
#include "DiscreteFissionPKAPDF.h"

template<>
InputParameters validParams<PKAGeneratorNeutronicsBase>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addClassDescription("PKA generator (neutronics) user object base class.\n Takes pdf and samples PKAs due to various interactions.");
  return params;
}

PKAGeneratorNeutronicsBase::PKAGeneratorNeutronicsBase(const InputParameters & parameters) :
    PKAGeneratorBase(parameters)
{
}
