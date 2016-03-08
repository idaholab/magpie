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
  Real nsq, x1, x2;

  // Marsaglia's method for uniformly sampling the surface of the sphere
  do
  {
    x1 = 2 * getRandomReal() - 1.0;
    x2 = 2 * getRandomReal() - 1.0;
    nsq = x1 * x1 + x2 * x2;
  } while (nsq >= 1);

  // construct normalized direction vector
  ion._dir(0) = 2.0 * x1 * std::sqrt(1.0 - nsq);
  ion._dir(1) = 2.0 * x2 * std::sqrt(1.0 - nsq);
  ion._dir(2) = 1.0 - 2.0 * nsq;
}
