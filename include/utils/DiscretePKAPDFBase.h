#ifndef DiscretePKAPDFBaseBASE_H
#define DiscretePKAPDFBaseBASE_H

#include "MultiIndex.h"
#include "mytrim/ion.h"

/**
 * Implements a discrete PDF for sampling
 * PKAs
 */
class DiscretePKAPDFBase
{
public:
  /// default constructor
  DiscretePKAPDFBase() {}

  /// constructor setting all necessary values
  DiscretePKAPDFBase(const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies);

  /// Uses the discrete probabilities for sampling the initial pka state
  virtual void drawSample(std::vector<MyTRIM_NS::IonBase> & initial_state) const = 0;

  /// accessor for getting the magnitude of this DiscretePKAPDF
  Real getMagnitude() const { return _magnitude; }

protected:
  /// NOTE: we pass by value here because we modify probabilities in the function for convenience
  virtual void precomputeCDF(MultiIndex<Real> probabilities) = 0;

  /// this method computes the magnitude but we cannot implement this in the base class
  virtual void computeMagnitude(MultiIndex<Real> probabilities) = 0;

  ///@{ helper function to draw a sample from a marginal probability function, returns "bin" id
  unsigned int sampleHelper(const MultiIndex<Real> & marginal_pdf) const;
  unsigned int sampleHelper(const MultiIndex<Real> & marginal_pdf, const std::vector<unsigned int> indices) const;
  unsigned int sampleHelper(const std::vector<Real> & marginal_pdf) const;
  ///@}

  /// magnitude for correct scaling with potential other DiscretePKAPDFBase objects
  Real _magnitude;

  /// vector storing the Z,A values in ZAID form
  std::vector<unsigned int> _zaids;

  /// number of zaids
  unsigned int _nZA;

  /// the energy group boundaries
  std::vector<Real> _energies;

  /// number of energy groups
  unsigned int _ng;
};

#endif // DiscretePKAPDFBase
