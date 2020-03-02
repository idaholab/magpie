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
 * Bulk Helmholtz energy psi^l, eq. (1) from Physical Review B 89, 184102 (2014)
 */
class PolarPFMPsiL : public DerivativeParsedMaterialHelper, public ExpressionBuilder
{
public:
  static InputParameters validParams();

  PolarPFMPsiL(const InputParameters & parameters);

protected:
  ///@{ Order parameter
  EBTerm _upsilon;
  EBTerm _theta;
  ///@}

  /// thermal energy of the melt
  const Real _g0;

  ///@{ interpolation coefficients
  const Real _a_theta;
  const Real _a_A;
  ///@}

  ///@{ change in thermal energy between phases
  const Real _deltaG10;
  const Real _deltaG21;
  ///@}

  ///@{ barrier coefficients
  const Real _a10;
  const Real _a20;
  const Real _a21;
  ///@}
};
