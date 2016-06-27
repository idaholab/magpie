#ifndef DiscretePKAPDFBaseBASE_H
#define DiscretePKAPDFBaseBASE_H

#include "MultiIndex.h"

/**
 * Implements a discrete PDF for sampling
 * PKAs
 */
class DiscretePKAPDFBase
{
public:
  DiscretePKAPDFBase(Real magnitude, const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies);

  /**
   * A struct storing the inital state of a primary knock-on atom
   */
  struct initialPKAState
  {
    unsigned int _Z;
    Real _mass;
    Real _energy;
    Point _direction;
  };

  /// Uses the discrete probabilities for sampling the initial pka state
  virtual void drawSample(std::vector<initialPKAState> & initial_state) = 0;

  /// NOTE: we pass by value here because we modify probabilities in the function for
  /// convenience
  virtual void precomputeCDF(MultiIndex<Real> probabilities) = 0;

protected:
  /// magnitude for correct scaling with potential other DiscretePKAPDFBase objects
  Real _magnitude;

  /// vector storing the Z,A values in ZAID form
  const std::vector<unsigned int> _zaids;
  /// number of zaids
  unsigned int _nZA;

  /// the energy group boundaries
  const std::vector<Real> _energies;
  /// number of energy groups
  unsigned int _ng;

  ///@{  helper function to draw a sample from a marginal probability function, returns "bin" id
  unsigned int sampleHelper(const MultiIndex<Real> & marginal_pdf) const;
  unsigned int sampleHelper(const MultiIndex<Real> & marginal_pdf, const std::vector<unsigned int> indices) const;
  unsigned int sampleHelper(const std::vector<Real> & marginal_pdf) const;
  ///@}

  /// uniformly samples direction
  std::vector<Real> sampleUniformDirection();
};

#endif // DiscretePKAPDFBase
