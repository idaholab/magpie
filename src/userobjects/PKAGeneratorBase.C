#include "PKAGeneratorBase.h"
#include "MagpieUtils.h"

template<>
InputParameters validParams<PKAGeneratorBase>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  return params;
}

PKAGeneratorBase::PKAGeneratorBase(const InputParameters & parameters) :
    DiscreteElementUserObject(parameters)
{
  setRandomResetFrequency(EXEC_TIMESTEP_END);
}

void
PKAGeneratorBase::setPosition(MyTRIM_NS::IonBase & ion) const
{
  ion._pos = MagpieUtils::randomElementPoint(*_current_elem, getRandomPoint());
}

void
PKAGeneratorBase::setRandomDirection(MyTRIM_NS::IonBase & ion) const
{
  const Real cosTheta = 2.0 * getRandomReal() - 1.0;
  const Real sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
  const Real phi = 2.0 * libMesh::pi * getRandomReal();
  ion._dir = Point(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}
