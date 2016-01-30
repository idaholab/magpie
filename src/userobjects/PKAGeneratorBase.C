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
PKAGeneratorBase::setPosition(MyTRIM_NS::ionBase & ion) const
{
  Point pos = MagpieUtils::randomElementPoint(*_current_elem, getRandomPoint());
  ion.pos[0] = pos(0);
  ion.pos[1] = pos(1);
  ion.pos[2] = pos(2);
}

void
PKAGeneratorBase::setRandomDirection(MyTRIM_NS::ionBase & ion) const
{
  Point dir;
  Real size_sq;

  do
  {
    dir = getRandomPoint() - Point(0.5, 0.5, 0.5);
    size_sq = dir.size_sq(); // soon norm_sq() ?

    // we reject points outside or the sphere with radius 1/2 (otherwise
    // there'd be a higehr probability to point towards the cube corners) and
    // points with small norms (for numerical reasons).
  } while (size_sq < 0.001 || size_sq > 0.25);

  // normalize direction vector
  dir /= std::sqrt(size_sq);

  ion.dir[0] = dir(0);
  ion.dir[1] = dir(1);
  ion.dir[2] = dir(2);
}
