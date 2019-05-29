/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#ifndef SCATTERINGRECOILCROSSSECTION_H
#define SCATTERINGRECOILCROSSSECTION_H

#include "GeneralUserObject.h"

class ScatteringRecoilCrossSection;

template <>
InputParameters validParams<ScatteringRecoilCrossSection>();

class ScatteringRecoilCrossSection : public GeneralUserObject
{
public:
  ScatteringRecoilCrossSection(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void finalize() override;

  /// number of neutron energy groups
  unsigned int nNeutronGroups() const { return _G; }

  /// number of recoil energy groups
  unsigned int nRecoilGroups() const { return _T; }

  /// the order of the legendre series expansion of the recoil cross section in mu_L
  unsigned int legendreOrder() const { return _L; }

  ///@{ neutron and recoil energy limits
  const std::vector<Real> & getNeutronEnergyLimits() const { return _neutron_energy_limits; };
  const std::vector<Real> & getRecoilEnergyLimits() const { return _recoil_energy_limits; };
  ///@}

  /// returns the recoil cross section Legendre expansion coefficient of order l for group g->t
  Real getSigmaRecoil(unsigned int g, unsigned int t, unsigned int l) const;

  ///@{ maximum and minimum laboratory scattering cosine for g->t group combination
  Real getMaxRecoilCosine(unsigned int g, unsigned t) const;
  Real getMinRecoilCosine(unsigned int g, unsigned t) const;
  ///@}

  /// is this g-> combination possible
  bool isRecoilPossible(unsigned int g, unsigned int t) const;

protected:
  /// helper function to get mu_C from T and E
  virtual Real getCMCosine(Real E, Real T, Real Q = 0.0) const = 0;

  /// helper function to get mu_L from mu_C & neutron energy
  virtual Real getLabCosine(Real E, Real T, Real Q = 0.0) const = 0;

  /// Method that finds the neutron energy group given a neutron energy
  unsigned int findNeutronEnergyGroup(Real energy);

  /// A wrapper method checking values printed to csv files
  Real csvPrint(Real value) const;

  /// below this number csv output is set to 0
  const Real _csv_tolerance;

  /// order of the quadrature
  unsigned int _quad_order;

  /// Function representing the neutron spectrum
  const Function & _neutron_spectrum;

  /// Order of Legendre polynomials
  unsigned int _L;

  /// Number of recoil energy bins
  unsigned int _T;

  /// Number of neutron energy groups
  unsigned int _G;

  /// Atomic mass of the isotope being analyzed
  const Real _atomic_mass;

  /// Mass fraction that limits the energy transfer from the neutron to the recoil
  Real _gamma;

  /// Mass fraction
  Real _alpha;

  /// Limits that characterize the neutron energy groups, must be in descending order
  const std::vector<Real> _neutron_energy_limits;

  /// Limits that characterize the recoil energy bins, must be in descending order
  const std::vector<Real> _recoil_energy_limits;

  /// Function representing the elastic cross section
  std::vector<const Function *> _scattering_cross_section;

  /// Varible for the neutron spectrum
  std::vector<Real> _xi_g;

  /// Elastic recoil cross section coefficients for the expansion in Legendre polynomials
  std::vector<std::vector<std::vector<Real>>> _recoil_coefficient;

  /// Varible to store the maximum and minimum possible cosines in the Lab frame
  std::vector<std::vector<std::vector<Real>>> _mu_L;

  /// Output file with the elastic recoil cross section coefficients for the expansion in Legendre polynomials
  std::string _recoil_xs_file_name;

  /// Output file with the maximum and minimum possible cosines in the Lab frame
  std::string _mu_L_file_name;

  ///@{ Gaussian quadrature rule
  std::vector<Real> _quad_points;
  std::vector<Real> _quad_weights;
  ///@}
};

#endif // SCATTERINGRECOILCROSSSECTION_H
#endif // GSL_ENABLED
