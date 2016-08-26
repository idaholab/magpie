#ifdef RATTLESNAKE_ENABLED
#ifndef NEUTRONICSSPECTRUMSAMPLERSN_H
#define NEUTRONICSSPECTRUMSAMPLERSN_H

#include "NeutronicsSpectrumSamplerBase.h"
#include "AQData.h"

// Forward Declarations
class NeutronicsSpectrumSamplerSN;

template<>
InputParameters validParams<NeutronicsSpectrumSamplerSN>();

/**
 * Computes PDFs from neutronics data that is used to sample PKAs due to
 * other reactions (not fission) for coupled BCMC simulations.
 * The reaction creating the PKAs _must_
 * match target_isotope_names with a unique recoil_isotope_names
 */
class NeutronicsSpectrumSamplerSN : public NeutronicsSpectrumSamplerBase
{
public:
  NeutronicsSpectrumSamplerSN(const InputParameters & parameters);

protected:
  /// a callback executed right before computeRadiatonDamagePDF
  virtual void preComputeRadiationDamagePDF();
  /// computes the PDF for isotope i, group g, and SH indices p
  virtual Real computeRadiationDamagePDF(unsigned int i, unsigned int g, unsigned int p);

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

#endif //NEUTRONICSSPECTRUMSAMPLERSN_H
#endif //RATTLESNAKE_ENABLED
