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
  Real norm_sq;

  do
  {
    ion._dir = getRandomPoint() - Point(0.5, 0.5, 0.5);
    norm_sq = ion._dir.norm_sq();

    // we reject points outside or the sphere with radius 1/2 (otherwise
    // there'd be a higher probability to point towards the cube corners) and
    // points with small norms (for numerical reasons).
  } while (norm_sq < 0.001 || norm_sq > 0.25);

  // normalize direction vector
  ion._dir /= std::sqrt(norm_sq);
}
