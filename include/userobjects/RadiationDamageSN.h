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
#ifdef RATTLESNAKE_ENABLED
#ifndef RADIATIONDAMAGESN_H
#define RADIATIONDAMAGESN_H

#include "RadiationDamageBase.h"

// Forward Declarations
class RadiationDamageSN;

template<>
InputParameters validParams<RadiationDamageSN>();

/**
 * Computes the PKA species/energy/direction distribution
 * at a given set of point for an SN calculation.
 * The reaction creating the PKAs _must_
 * match target_isotope_names with a unique recoil_isotope_names
 */
class RadiationDamageSN : public RadiationDamageBase
{
public:
  RadiationDamageSN(const InputParameters & parameters);

protected:
  /// a callback executed right before computePKA
  virtual void preComputePKA();
  /// computes the PKA for isotope i, group g, and SH indices p
  virtual Real computePKA(unsigned int i, unsigned int g, unsigned int p);

  /// vector of target zaids
  const std::vector<std::string> & _recoil_isotope_names;
  /// angular quadrature object
  const AngularQuadrature & _aq;
  /// Number of angular directions in _aq
  unsigned int _ndir;
  /// spherical harmonics coefficient object
  const SHCoefficients _shm;
  /// the angular flux
  std::vector<std::vector<const VariableValue *> > _angular_flux;
  /// the angular flux moment evaluated at one quadrature point
  std::vector<std::vector<Real> > _flux_moment;
  /// stores the recoil cross sections
  std::vector<std::vector<std::vector<std::vector<std::vector<Real> > > > > _recoil_cross_section;
};

#endif //RADIATIONDAMAGESN_H
#endif //RATTLESNAKE_ENABLED
