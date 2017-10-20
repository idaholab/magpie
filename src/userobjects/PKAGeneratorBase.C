/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAGeneratorBase.h"
#include "MagpieUtils.h"
#include <algorithm>

template<>
InputParameters validParams<PKAGeneratorBase>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  params.addClassDescription("PKA generator user object base class.\n Takes pdf and samples PKAs due to various interactions.");
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

int
PKAGeneratorBase::ionTag(const std::vector<Real> & rasterizer_Z, const std::vector<Real> & rasterizer_m, Real Z, Real m) const
{
  // this function relies on the exact representation of whole numbers in IEEE floating point numbers
  // up to a reasonable upper limit [Z < m < 300]
  const auto & it = std::find(rasterizer_Z.begin(), rasterizer_Z.end(), Z);
  if (it != rasterizer_Z.end())
  {
    unsigned int index = std::distance(rasterizer_Z.begin(), it);
    if (rasterizer_m[index] == m)
      return index;
  }
  return -1;
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
