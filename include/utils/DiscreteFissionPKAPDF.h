#ifndef DISCRETEFISSIONPKAPDF_H
#define DISCRETEFISSIONPKAPDF_H

#include "DiscretePKAPDFBase.h"
#include "MultiIndex.h"
#include "MooseEnum.h"
#include "MagpieUtils.h"

/**
 * Implements a discrete PDF for sampling fission products given target, energy,
 * and fission rates per target
 */
class DiscreteFissionPKAPDF : DiscretePKAPDFBase
{
public:
  DiscreteFissionPKAPDF(Real magnitude,const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities);

  /// override drawSample
  virtual void drawSample(std::vector<initialPKAState> & initial_state) override;

  /// override preComputeCDF. NOTE: we pass by value here because we modify probabilities in the function for
  /// convenience
  virtual void precomputeCDF(MultiIndex<Real> probabilities) override;

  /// reads fission yield data
  void readFissionData(const std::vector<unsigned int> & ZAID);

protected:
  /**
   * reads Z and A number and returns kinetic energy in MeV
   */
  Real determineFragmentsEnergy(unsigned int Z, unsigned int A);

  /// samples the number neutrons per fission based on target and energy
  unsigned int sampleNu(MagpieUtils::NeutronEnergyType energy_type, unsigned int zaid);

  /// we can cache the marginal distributions
  MultiIndex<Real> _marginal_cdf_target;
  MultiIndex<Real> _conditional_cdf_energy;

  /**
   * vector of maps of fission product zaids and cdf
   * energy type enum -> map (zaid -> vector)
   */
  std::vector<std::map<unsigned int, std::vector<unsigned int> > >  _fission_zaids;
  std::vector<std::map<unsigned int, std::vector<Real> > >  _fission_cdf;
};

#endif //DISCRETEFISSIONPKAPDF_H
