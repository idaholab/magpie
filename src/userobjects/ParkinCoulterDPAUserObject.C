/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "ParkinCoulterDPAUserObject.h"
#include "PolyatomicRecoil.h"
#include "PolyatomicDisplacementFunction.h"
#include "PolyatomicDamageEnergyFunction.h"
#include "PolyatomicDisplacementDerivativeFunction.h"
#include "MooseMesh.h"

// various other includes
#include <mytrim/element.h>
#include <limits>

registerMooseObject("MagpieApp", ParkinCoulterDPAUserObject);

InputParameters
ParkinCoulterDPAUserObject::validParams()
{
  InputParameters params = DPAUserObjectBase::validParams();
  params += ParkinCoulterInterface::validParams();
  params.addClassDescription(
      "Computes the dose in dpa from composition, cross section, damage type, and neutron flux for "
      "polyatomic materials using Parkin-Coulter's method.");
  return params;
}

ParkinCoulterDPAUserObject::ParkinCoulterDPAUserObject(const InputParameters & parameters)
  : DPAUserObjectBase(parameters), ParkinCoulterInterface(this)
{
}

void
ParkinCoulterDPAUserObject::initialSetup()
{
  prepare();
}

std::vector<unsigned int>
ParkinCoulterDPAUserObject::atomicNumbers() const
{
  return getAtomicNumbers();
}

std::vector<Real>
ParkinCoulterDPAUserObject::massNumbers() const
{
  return getMassNumbers();
}

std::vector<Real>
ParkinCoulterDPAUserObject::numberFractions() const
{
  return getNumberFractions();
}

void
ParkinCoulterDPAUserObject::initDamageFunctions()
{
  _padf = libmesh_make_unique<PolyatomicDisplacementFunction>(polyMat(), NET, _Ecap);
}

void
ParkinCoulterDPAUserObject::execute()
{
  accumulateDamage();
}

void
ParkinCoulterDPAUserObject::onCompositionChanged()
{
  computeDamageFunctions();
  _padf->computeDisplacementFunctionIntegral();
}

Real
ParkinCoulterDPAUserObject::integralDamageFunction(Real T, unsigned int i, unsigned int j) const
{
  return _padf->linearInterpolationIntegralDamageFunction(T, i, j);
}

void
ParkinCoulterDPAUserObject::finalize()
{
}

#endif
