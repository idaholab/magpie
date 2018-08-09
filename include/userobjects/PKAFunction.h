/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef PKAFUNCTION_H
#define PKAFUNCTION_H

#include "PKAEmpiricalBase.h"
#include "Function.h"

class PKAFunction;

template <>
InputParameters validParams<PKAFunction>();

/**
 * PKAs with time dependent mass, charge, energy, and rate
 */
class PKAFunction : public PKAEmpiricalBase
{
public:
  PKAFunction(const InputParameters & parameters);

protected:
  /// Fission rate (per unit volume)
  virtual Real getPKARate() const override { return _pka_rate.value(_time, Point()); };

  ///@{ charge, mass, energy
  virtual unsigned int getZ() const override { return _Z.value(_time, Point()); };
  virtual Real getM() const override { return _m.value(_time, Point()); };
  virtual Real getE() const override { return _E.value(_time, Point()); };
  ///@}

  /// Fission rate (per unit volume)
  Function & _pka_rate;

  /// PKA nuclear charge
  Function & _Z;

  /// PKA mass
  Function & _m;

  /// PKA Energy (in eV)
  Function & _E;

  /// time
  const Real & _time;
};

#endif // PKAFUNCTION_H
