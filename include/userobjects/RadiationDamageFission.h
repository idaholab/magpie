/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#if RATTLESNAKE_ENABLED
#ifndef RADIATIONDAMAGEFISSION_H
#define RADIATIONDAMAGEFISSION_H

#include "RadiationDamageBase.h"

// Forward Declarations
class RadiationDamageFission;

template<>
InputParameters validParams<RadiationDamageFission>();

/**
 * Computes the PKA species/energy/direction distribution
 * for fission reactions.
 * NOTE: Currently fission is assumed to be isotropic in the LAB
 * frame regardless of the incoming energy.
 */
class RadiationDamageFission : public RadiationDamageBase
{
public:
  RadiationDamageFission(const InputParameters & parameters);

protected:
  /// computes the PKA for isotope i, group g, and SH indices p
  /// NOTE: for fission p does not mateter
  virtual Real computePKA(unsigned int i, unsigned int g, unsigned int /*p*/);

  /// the angular flux
  std::vector<const VariableValue *> _scalar_flux;
  /// stores the recoil cross sections
  std::vector<std::vector<std::vector<Real> > > _fission_cross_section;
};

#endif //RADIATIONDAMAGEFISSION_H
#endif //RATTLESNAKE_ENABLED
