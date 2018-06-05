/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKAGENERATORNEUTRONICSBASE_H
#define PKAGENERATORNEUTRONICSBASE_H

#include "PKAGeneratorBase.h"
#include "MultiIndex.h"
#include "MyTRIMRasterizer.h"
#include "mytrim/ion.h"
#include "DiscreteFissionPKAPDF.h"

class PKAGeneratorNeutronicsBase;

template<>
InputParameters validParams<PKAGeneratorNeutronicsBase>();

/**
 * A PKAGenerator class that uses neutronics computed reaction rates for
 * sampling the inital state of PKAs.
 */
class PKAGeneratorNeutronicsBase : public PKAGeneratorBase
{
public:
  PKAGeneratorNeutronicsBase(const InputParameters & parameters);

  /// helper function to set pdf based on neutronics information from Mammoth
  virtual void setPDF(const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities) = 0;

protected:
  /**
   * The partial reaction rates are reaction rates of a certain type from the neutronics calculation
   * for a single nculide, e.g. fission rate of U-235. NOTE: it multiplies the smoothly varying _average_
   * of the U-235 number density already
   *
   * In Magpie, the rate of PKA creation is obtained by multiplying with the _local_ number density. To obtain a
   * proper PKA creation rate, we need to divide by the smooth average. This is what _averaged_number_densities is for.
   */
  std::vector<const Real *> _partial_neutronics_reaction_rates;
  std::vector<Real> _stored_reaction_rates;
  std::vector<const Real *> _averaged_number_densities;
  std::vector<Real> _stored_densities;
};

#endif // PKAGENERATORNEUTRONICSBASE_H
