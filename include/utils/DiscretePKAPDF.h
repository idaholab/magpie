#ifndef DISCRETEPKAPDF_H
#define DISCRETEPKAPDF_H

#include "DiscretePKAPDFBase.h"
#include "MultiIndex.h"

/**
 * Implements a discrete PDF for sampling
 * PKAs
 */
class DiscretePKAPDF : DiscretePKAPDFBase
{
public:
  DiscretePKAPDF(Real magnitude, std::vector<unsigned int> ZAID, std::vector<Real> energies,
                 unsigned int na, unsigned int np, MultiIndex<Real> _probabilities);

  /// override drawSample
  void drawSample(initialPKAState & initial_state);

protected:
  /// helper functions that precomputes marginal cdfs
  void precomputeCDF();

  ///@{  helper function to draw a sample from a marginal probability function, returns "bin" id
  unsigned int sampleHelper(const MultiIndex<Real> & marginal_pdf) const;
  unsigned int sampleHelper(const MultiIndex<Real> & marginal_pdf, const std::vector<unsigned int> indices) const;
  ///@}

  /// vector storing the Z,A values in ZAID form
  const std::vector<unsigned int> _zaids;
  /// number of zaids
  unsigned int _nZA;

  /// the energy group boundaries
  const std::vector<Real> _energies;
  /// number of energy groups
  unsigned int _ng;

  /// number of azimuthal bins
  unsigned int _na;
  /// vector storing the azimuthal angle boundaries
  const Real _dphi;

  /// number of polar bins
  unsigned int _np;
  /// polar cosine bin width
  const Real _dmu;

  /// stores the discrete probabilities for (Z,A)/Energy/Direction(phi, mu)
  MultiIndex<Real> _probabilities;

  /// we can cache the marginal distributions
  MultiIndex<Real> _marginal_cdf_mu;
  MultiIndex<Real> _marginal_cdf_phi;
  MultiIndex<Real> _marginal_cdf_energy;
  MultiIndex<Real> _marginal_cdf_zaid;
};

#endif // DiscretePKAPDF
