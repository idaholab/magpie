/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKAFISSIONFRAGMENTNEUTRONICS_H
#define PKAFISSIONFRAGMENTNEUTRONICS_H

#include "PKAGeneratorBase.h"
#include "PKAGeneratorNeutronicsBase.h"
#include "DiscreteFissionPKAPDF.h"

class PKAFissionFragmentNeutronics;

template<>
InputParameters validParams<PKAFissionFragmentNeutronics>();

class PKAFissionFragmentNeutronics : public PKAGeneratorNeutronicsBase
{
public:
  PKAFissionFragmentNeutronics(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData &) const;

  virtual void setPDF(const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities);

protected:
  /**
   * Given the magnitude, ZAIDs, energies, and probabilities passed from MAMMOTH,
   * the cdf will be precomputed for each fission fragment and energy. These cdfs
   * will be sampled to determine the initial state of the ions passed to MyTRIM
   */
  DiscreteFissionPKAPDF _pdf;

  /**
   * _partial_fission_rates: the fission rate per nuclide is:
   * f_i = N_i * (sum_g sigma_{f,g,i} phi_g),
   * where i is the nuclide ID and g is the energy group.
   * N_i is provided as variable in the Rasterizer, but the (sum_g sigma_{f,g,i} phi_g) must
   * be provided as partial fission rate here
   * Note: the fission rate densities are in units of [ fissions / ( mesh-length-unit^3 sec ) ]
   */
  std::vector<const Real *> _partial_fission_rates;
  std::vector<Real> _stored_pps;
};

#endif //PKAFISSIONFRAGMENTNEUTRONICS_H
