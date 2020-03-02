/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "PKAEmpiricalBase.h"

/**
 * PKAs with constant mass, charge, energy, and rate
 */
class PKAConstant : public PKAEmpiricalBase
{
public:
  static InputParameters validParams();

  PKAConstant(const InputParameters & parameters);

protected:
  /// Fission rate (per unit volume)
  virtual Real getPKARate() const override { return _pka_rate; };

  ///@{ charge, mass, energy
  virtual unsigned int getZ() const override { return _Z; };
  virtual Real getM() const override { return _m; };
  virtual Real getE() const override { return _E; };
  ///@}

  /// Fission rate (per unit volume)
  const Real _pka_rate;

  /// PKA nuclear charge
  const unsigned int _Z;

  /// PKA mass
  const Real _m;

  /// PKA Energy (in eV)
  const Real _E;
};
