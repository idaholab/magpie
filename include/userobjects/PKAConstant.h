/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKACONSTANT_H
#define PKACONSTANT_H

#include "PKAGeneratorBase.h"

class PKAConstant;

template<>
InputParameters validParams<PKAConstant>();

/**
 *
 */
class PKAConstant : public PKAGeneratorBase
{
public:
  PKAConstant(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData & averaged_data) const;

protected:
  /// Fission rate (per unit volume)
  const Real _pka_rate;

  /// PKA nuclear charge
  const unsigned int _Z;

  /// PKA mass
  const Real _m;

  /// PKA Energy (in eV)
  const Real _E;
};

#endif // PKACONSTANT_H
