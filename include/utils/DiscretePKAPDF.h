/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

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
  /// default constructor
  DiscretePKAPDF();

  DiscretePKAPDF(const std::vector<unsigned int> & ZAID,
                 const std::vector<Real> & energies,
                 const MultiIndex<Real> & probabilities);

  /// override drawSample
  virtual void drawSample(std::vector<MyTRIM_NS::IonBase> & initial_state) const override;

  /// overload the outstream operator
  friend std::ostream & operator<<(std::ostream & out, const DiscretePKAPDF & pdf);

protected:
  /// NOTE: we pass by value here because we modify probabilities in the function for convenience
  virtual void precomputeCDF(MultiIndex<Real> probabilities) override;

  /// this method computes the magnitude encoded in probabilities
  virtual void computeMagnitude(MultiIndex<Real> probabilities) override;

  /// number of azimuthal bins
  unsigned int _na;

  /// vector storing the azimuthal angle boundaries
  Real _dphi;

  /// number of polar bins
  unsigned int _np;

  /// polar cosine bin width
  Real _dmu;

  /// store the pdf
  MultiIndex<Real> _probability_density_function;

  /// we can cache the marginal distributions
  MultiIndex<Real> _marginal_cdf_mu;
  MultiIndex<Real> _marginal_cdf_phi;
  MultiIndex<Real> _marginal_cdf_energy;
  MultiIndex<Real> _marginal_cdf_zaid;
};

#endif // DiscretePKAPDF
