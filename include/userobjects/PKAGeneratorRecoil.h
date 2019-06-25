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
#include "DiscretePKAPDF.h"

class PKAGeneratorRecoil;

template <>
InputParameters validParams<PKAGeneratorRecoil>();

class PKAGeneratorRecoil : public PKAGeneratorNeutronicsBase
{
public:
  PKAGeneratorRecoil(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> &,
                          const MyTRIMRasterizer::PKAParameters &,
                          const MyTRIMRasterizer::AveragedData &) const;

  virtual void setPDF(const std::vector<unsigned int> & ZAID,
                      const std::vector<Real> & energies,
                      const MultiIndex<Real> & probabilities);

protected:
  /**
   * Discrete PDF object from which the nuclide, energy, and direction of recoils
   * are sampled
   */
  DiscretePKAPDF _pdf;
};

#endif // PKAGeneratorRecoil_H
