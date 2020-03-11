/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include "GeneralUserObject.h"

class DPAUserObjectBase : public GeneralUserObject
{
public:
  static InputParameters validParams();

  DPAUserObjectBase(const InputParameters & parameters);

  ///@{ get the
  Real getDPA() const { return _dpa; }
  Real getPartialDPA(unsigned int Z, Real A) const;
  ///@}

  void initialize() override {}

protected:
  void prepare();
  bool changed() const;
  std::vector<unsigned int> getAtomicNumbers() const;
  std::vector<Real> getMassNumbers() const;
  std::vector<Real> getNumberFractions() const;
  Real getMaxEnergy() const;

  /// accumulated dpa
  void accumulateDamage();

  /// a helper function that sets _ns and checks consistency of A, Z, N
  void initAZNHelper();

  /// a helper function that computes the neutron damage efficiency
  Real
  neutronDamageEfficiency(unsigned int i, unsigned int j, unsigned int g, unsigned int x) const;

  /// a virtual function computing the integral damage function
  virtual Real integralDamageFunction(Real T, unsigned int i, unsigned int j) const = 0;

  /// callback that is executed when composition changed and damage functions must be recomputed
  virtual void onCompositionChanged() = 0;

  /// a helper that assigns a unique string to a Z, A pair
  std::string zaidHelper(unsigned int Z, Real A) const;

  /// tolerance for recomputing the displacement function
  Real _tol = 1e-10;

  /// the computed dose rates
  Real _dpa = 0;

  /// the computed dose rate by species; this is a map because presence of (Z,A) can change dynamically
  std::map<std::string, Real> _partial_dpa;

  /// is damage accumulated during a transient or computed for steady state
  bool _is_transient_irradiation;

  /// irradiation_time used when dpa is estimated from steady-state calculations
  Real _irradiation_time;

  /// the neutron reaction types considered for computing dpa
  MultiMooseEnum _neutron_reaction_types;

  /// number of reaction types creating radiation damage
  unsigned int _nr;

  ///@{ data used for computing dpa value
  std::vector<Real> _atomic_numbers;
  std::vector<Real> _mass_numbers;
  std::vector<Real> _number_densities;
  std::vector<Real> _energy_group_boundaries;
  std::vector<Real> _scalar_flux;
  std::vector<std::vector<std::vector<Real>>> _cross_sections;
  ///@}

  /// Q values for each reaction type and isotope
  std::vector<std::vector<Real>> _q_values;

  ///@{ the "old" versions of the data; used for determining if disp function update is required
  std::vector<Real> _atomic_numbers_old;
  std::vector<Real> _mass_numbers_old;
  std::vector<Real> _number_densities_old;
  ///@}

  /// number of neutron energy groups
  unsigned int _ng;
};

#endif // GSL_ENABLED
