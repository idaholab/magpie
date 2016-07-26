#ifndef DISCRETEPKAPDF_H
#define DISCRETEPKAPDF_H

#include "DiscretePKAPDFBase.h"
#include "MultiIndex.h"
#include "mytrim/ion.h"

/**
 * Implements a discrete PDF for sampling
 * PKAs
 */
class DiscretePKAPDF : DiscretePKAPDFBase
{
public:
  DiscretePKAPDF(Real magnitude, const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities);

  /// override drawSample
  virtual void drawSample(std::vector<MyTRIM_NS::IonBase> & initial_state) const override;

  /// NOTE: we pass by value here because we modify probabilities in the function for
  /// convenience
  virtual void precomputeCDF(MultiIndex<Real> probabilities) override;

protected:


  /// number of azimuthal bins
  unsigned int _na;
  /// vector storing the azimuthal angle boundaries
  const Real _dphi;

  /// number of polar bins
  unsigned int _np;
  /// polar cosine bin width
  const Real _dmu;

  /// we can cache the marginal distributions
  MultiIndex<Real> _marginal_cdf_mu;
  MultiIndex<Real> _marginal_cdf_phi;
  MultiIndex<Real> _marginal_cdf_energy;
  MultiIndex<Real> _marginal_cdf_zaid;
};

#endif // DiscretePKAPDF
