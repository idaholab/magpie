/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "DiscretePKAPDFBase.h"
#include "MultiIndex.h"
#include "MagpieUtils.h"
#include "mytrim/ion.h"

/**
 * Implements a discrete PDF for sampling fission products given target, energy,
 * and fission rates per target
 */
class DiscreteFissionPKAPDF : DiscretePKAPDFBase
{
public:
  /// the default constructor
  DiscreteFissionPKAPDF();

  /// the actual constructor
  DiscreteFissionPKAPDF(const std::vector<unsigned int> & ZAID,
                        const std::vector<Real> & energies,
                        const MultiIndex<Real> & probabilities);

  /// override drawSample
  virtual void drawSample(std::vector<MyTRIM_NS::IonBase> & initial_state) const override;

protected:
  /// override preComputeCDF. NOTE: we pass by value here because we modify probabilities in the function for convenience
  virtual void precomputeCDF(MultiIndex<Real> probabilities) override;

  /// reads fission yield data
  void readFissionData(const std::vector<unsigned int> & ZAID);

  /// this method computes the magnitude encoded in probabilities
  virtual void computeMagnitude(MultiIndex<Real> probabilities) override;

  /**
   * reads Z and A number and returns kinetic energy in MeV
   */
  Real determineFragmentsEnergy(unsigned int Z, unsigned int A) const;

  /// samples the number neutrons per fission based on target and energy
  unsigned int sampleNu(MagpieUtils::NeutronEnergyType energy_type, unsigned int zaid) const;

  /// we can cache the marginal distributions
  MultiIndex<Real> _marginal_cdf_target;
  MultiIndex<Real> _conditional_cdf_energy;

  /**
   * vector of maps of fission product zaids and cdf
   * energy type enum -> map (zaid -> vector)
   */
  std::vector<std::map<unsigned int, std::vector<unsigned int>>> _fission_zaids;
  std::vector<std::map<unsigned int, std::vector<Real>>> _fission_cdf;
};
