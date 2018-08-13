/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKAGENERATORALPHADECAY_H
#define PKAGENERATORALPHADECAY_H

#include "PKAGeneratorBase.h"

class PKAGeneratorAlphaDecay;

template <>
InputParameters validParams<PKAGeneratorAlphaDecay>();

/**
 * PKAGenerator for alpha decays and corresponding recoil atoms
 */
class PKAGeneratorAlphaDecay : public PKAGeneratorBase
{
public:
  PKAGeneratorAlphaDecay(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> &,
                          const MyTRIMRasterizer::PKAParameters &,
                          const MyTRIMRasterizer::AveragedData &) const override;

  /// this stores the decay information for a single alpha decaying nuclide
  struct DecayData
  {
    Real _decay_constants;
    Real _alpha_energies;
    Real _intensities;
  };

protected:
  void readAlphaData();

  /// the Z/A ids required for this calculation
  const std::vector<unsigned int> _zaids;

  /// the decay data for each zaid
  std::map<unsigned int, std::vector<DecayData>> _decay_data_sets;

  /// half-life conversion from seconds to model unit
  Real _time_conversion;
};

#endif // PKAGENERATORALPHADECAY_H
