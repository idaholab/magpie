#include "PKAGeneratorBase.h"

template<>
InputParameters validParams<PKAGeneratorBase>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  return params;
}

PKAGeneratorBase::PKAGeneratorBase(const InputParameters & parameters) :
    DiscreteElementUserObject(parameters)
{
}
