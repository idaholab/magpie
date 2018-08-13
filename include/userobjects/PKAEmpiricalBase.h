/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKAEMPIRICALBASE_H
#define PKAEMPIRICALBASE_H

#include "PKAGeneratorBase.h"

class PKAEmpiricalBase;

template <>
InputParameters validParams<PKAEmpiricalBase>();

/**
 * Base class for empirical PKA generators
 */
class PKAEmpiricalBase : public PKAGeneratorBase
{
public:
  PKAEmpiricalBase(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> &,
                          const MyTRIMRasterizer::PKAParameters &,
                          const MyTRIMRasterizer::AveragedData &) const;

protected:
  /// Fission rate (per unit volume)
  virtual Real getPKARate() const = 0;

  ///@{ charge, mass, energy
  virtual unsigned int getZ() const = 0;
  virtual Real getM() const = 0;
  virtual Real getE() const = 0;
  ///@}
};

#endif // PKAEMPIRICALBASE_H
