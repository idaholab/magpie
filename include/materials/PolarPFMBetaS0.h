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
  EBTerm _theta;
  const Real _a_beta;
  const Real _beta10;
  const Real _beta20;
};
