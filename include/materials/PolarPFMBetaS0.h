/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "DerivativeParsedMaterialHelper.h"
#include "ExpressionBuilder.h"

/**
 * Solid-melt gradient energy coefficient, eq. (7) from Physical Review B 89, 184102 (2014)
 */
class PolarPFMBetaS0 : public DerivativeParsedMaterialHelper, public ExpressionBuilder
{
public:
  static InputParameters validParams();

  PolarPFMBetaS0(const InputParameters & parameters);

protected:
  /// Coupled solid phase order parameter
  EBTerm _theta;

  /// interpolation function parameter
  const Real _a_beta;

  /// solid1 / melt gradient energy coefficient
  const Real _beta10;
  /// solid2 / melt gradient energy coefficient
  const Real _beta20;
};
