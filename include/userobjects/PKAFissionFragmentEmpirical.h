/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKAFISSIONFRAGMENTEMPIRICAL_H
#define PKAFISSIONFRAGMENTEMPIRICAL_H

#include "PKAGeneratorBase.h"
#include "mytrim/invert.h"

class PKAFissionFragmentEmpirical;

template <>
InputParameters validParams<PKAFissionFragmentEmpirical>();

/**
 *
 */
class PKAFissionFragmentEmpirical : public PKAGeneratorBase
{
public:
  PKAFissionFragmentEmpirical(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                          Real dt,
                          Real vol,
                          Real recoil_rate_scaling,
                          const MyTRIMRasterizer::AveragedData & averaged_data) const;

protected:
  /// Fission rate (per unit volume) assuming pure fully dense UO2
  const PostprocessorValue & _fission_rate;

  /**
   * Variable for the relative Uranium density (0..1).
   * Use a CONSTANT MONOMIAL AuxVariable for this.
   */
  const VariableValue & _relative_density;
};

#endif // PKAFISSIONFRAGMENTEMPIRICAL_H
