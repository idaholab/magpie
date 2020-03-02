/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "PKAGeneratorBase.h"
#include "PKAGeneratorNeutronicsBase.h"
#include "DiscreteFissionPKAPDF.h"

class PKAFissionFragmentNeutronics : public PKAGeneratorNeutronicsBase
{
public:
  static InputParameters validParams();

  PKAFissionFragmentNeutronics(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> &,
                          const MyTRIMRasterizer::PKAParameters &,
                          const MyTRIMRasterizer::AveragedData &) const;

  virtual void setPDF(const std::vector<unsigned int> & ZAID,
                      const std::vector<Real> & energies,
                      const MultiIndex<Real> & probabilities);

protected:
  /**
   * Given the magnitude, ZAIDs, energies, and probabilities passed from MAMMOTH,
   * the cdf will be precomputed for each fission fragment and energy. These cdfs
   * will be sampled to determine the initial state of the ions passed to MyTRIM
   */
  DiscreteFissionPKAPDF _pdf;
};
