#include "PKAGeneratorNeutronicsBase.h"
#include "MagpieUtils.h"
#include "MultiIndex.h"
#include "DiscreteFissionPKAPDF.h"

template<>
InputParameters validParams<PKAGeneratorNeutronicsBase>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  return params;
}

PKAGeneratorNeutronicsBase::PKAGeneratorNeutronicsBase(const InputParameters & parameters) :
    PKAGeneratorBase(parameters)
{
}
