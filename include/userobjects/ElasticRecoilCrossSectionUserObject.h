
#ifndef ELASTICRECOILCROSSSECTIONUSEROBJECT_H
#define ELASTICRECOILCROSSSECTIONUSEROBJECT_H

#include "GeneralUserObject.h"

class ElasticRecoilCrossSectionUserObject;

template <>
InputParameters validParams<ElasticRecoilCrossSectionUserObject>();

class ElasticRecoilCrossSectionUserObject : public GeneralUserObject
{
public:
  ElasticRecoilCrossSectionUserObject(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override;

  virtual void finalize() override;

  ///@{ accessors for recoil cross section data
  unsigned int nNeutronGroups() const { return _G; }
  unsigned int nRecoilGroups() const { return _T; }
  unsigned int legendreOrder() const { return _L; }
  const std::vector<Real> & getNeutronEnergyLimits() const { return _neutron_energy_limits; };
  const std::vector<Real> & getRecoilEnergyLimits() const {return _recoil_energy_limits; };
  Real getSigmaRecoil(unsigned int g, unsigned int t, unsigned int l) const;
  Real getMaxRecoilCosine(unsigned int g, unsigned t) const;
  Real getMinRecoilCosine(unsigned int g, unsigned t) const;
  ///@}

protected:
  /// Method that finds the neutron energy group given a neutron energy
  unsigned int findNeutronEnergyGroup(Real energy);

  /// A wrapper method checking values printed to csv files
  Real csvPrint(Real value) const;

  /// below this number csv output is set to 0
  const Real _csv_tolerance;

  /// order of the quadrature
  unsigned int _quad_order;

  /// Function representing the neutron spectrum
  Function & _neutron_spectrum;

  /// Function representing the scattering law
  Function & _scattering_law;

  /// Function representing the elastic cross section
  Function & _elastic_xs;

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

  /// Varible for the neutron spectrum
  std::vector<Real> _xi_g;

  /// Elastic recoil cross section coefficients for the expansion in Legendre polynomials
  std::vector<std::vector<std::vector<Real>>> _erxs_coeff;

  /// Varible to store the maximum and minimum possible cosines in the Lab frame
  std::vector<std::vector<std::vector<Real>>> _save_mu_L;

  /// Output file with the elastic recoil cross section coefficients for the expansion in Legendre polynomials
  std::string _erxs_file_name;

  /// Output file with the maximum and minimum possible cosines in the Lab frame
  std::string _mu_L_file_name;

  ///@{ Gaussian quadrature rule
  std::vector<Real> _quad_points;
  std::vector<Real> _quad_weights;
  ///@}
};

#endif
