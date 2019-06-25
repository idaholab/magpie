/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <mytrim/material.h>

// forward declarations
class SimconfType;
class IonBase;

enum nrt_type
{
  ENERGY,
  TOTAL,
  NET,
  NET_DERIVATIVE
};

/**
 * Implements the computation of polyatomic displacement functions using Parkin & Coulter's
 * method
 */
class PolyatomicDisplacementFunctionBase
{
public:
  /// default constructor
  PolyatomicDisplacementFunctionBase(std::vector<MyTRIM_NS::Element> polyatomic_material,
                                     nrt_type damage_function_type,
                                     std::vector<std::vector<Real>> Ecap);

  virtual ~PolyatomicDisplacementFunctionBase();

  void advanceDisplacements(Real Emax);

  ///@{ some getters needed for accessing this pointer in odeRHS
  unsigned int nSpecies() const { return _n_species; }
  unsigned int problemSize() const { return _problem_size; }
  unsigned int nEnergySteps() const { return _energy_history.size(); }
  unsigned int findSpeciesIndex(unsigned int atomic_number, Real mass_number) const;
  Real minEnergy() const { return _energy_history[0]; }
  Real energyPoint(unsigned int na) const { return _energy_history[na]; }
  Real lambda(unsigned int i, unsigned int j) const { return _lambda[i][j]; }
  Real numberFraction(unsigned int i) const { return _material->_element[i]._t; }
  Real numberDensity(unsigned int i) const { return _material->_element[i]._t * _material->_arho; }
  ///@}

  /// gets stopping power for a given species and energy; non-const because it uses _ions so no need to construct ion
  Real stoppingPower(unsigned int species, Real energy);

  /// returns derivative of the stopping power of species w.r.t. to changes in the number fraction of changing_species
  Real stoppingPowerDerivative(unsigned int species, unsigned int changing_species, Real energy);

  /// retrives the grid index for an energy value
  unsigned int energyIndex(Real energy) const;

protected:
  /// computes the integral int_0^t dT T * d(sigma_ij) / dT for species combination i, j and small t
  Real
  weightedScatteringIntegral(Real energy, Real energy_limit, unsigned int i, unsigned int j) const;

  /// atomic scattering cross section d sigma_ij / dT, i: projectile, j: target
  Real
  scatteringCrossSection(unsigned int i, unsigned int j, Real energy, Real recoil_energy) const;

  /// the universal scattering function by Lindhard usually denoted as f(t^1/2)
  Real universalF(Real xi) const;

  /// displacement probability
  virtual Real displacementProbability(unsigned int k, Real energy) const;

  /// capture probability
  virtual Real
  nonCaptureProbability(unsigned int i, unsigned int k, Real energy, Real recoil_energy) const;

  /// damage function type [nij and gij, respectively in PK JNM 101, 1981; or nu_i JNM 88, (1980)]
  nrt_type _damage_function_type;

  /// order of the quadrature
  unsigned int _quad_order;

  /// the number of different species in the material
  unsigned int _n_species;

  /// the size of the problem = _n_species**2 for TOTAL & NET, _n_species for ENERGY
  unsigned int _problem_size;

  /// the current maximum energy to which PC equations are solved
  std::vector<Real> _energy_history;

  /**
   * The displacement values nu_ij flattened into 1D array
   */
  std::vector<std::vector<Real>> _displacement_function;

  /// Ecap_ik average residual energy of type i atom to not be trapped at k-site
  std::vector<std::vector<Real>> _Ecap;

  /// _lambda = 4.0 * Ai * Aj / (Ai + Aj) / (Ai + Aj)
  std::vector<std::vector<Real>> _lambda;

  /// Bohr radius in A
  const Real _abohr = 0.529177;

  /**
   * energy threshold below which the integral over cross sections is treated using
   * the asymptotic form of the cross sections for low energies and nu_j(T) is replaced
   * by nu_j(T) \approx T.
   */
  Real _asymptotic_threshold = 0.1;

  /// integral of the scattering cross sections int_{Edisp}^{Lambda_ij E} dT d(sigma_ij)/dT; Edisp: target displacement threshold
  // std::vector<std::vector<Real>> _scattering_xs_integral;

  ///@{ MYTRIM objects to use compute stopping power and cross sections
  std::unique_ptr<MyTRIM_NS::SimconfType> _simconf;
  std::unique_ptr<MyTRIM_NS::MaterialBase> _material;
  std::vector<std::unique_ptr<MyTRIM_NS::IonBase>> _ions;
  ///@}

  ///@{ GSL ODE system and driver to solve ODE
  gsl_odeiv2_system _sys;
  gsl_odeiv2_driver * _ode_driver;
  ///@}

  ///@{ Gaussian quadrature rule
  std::vector<Real> _quad_points;
  std::vector<Real> _quad_weights;
  ///@}
};

#endif // DiscretePKAPDF
#endif
