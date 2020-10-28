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
 * Interpolating function eq. (9) with a prefactor of beta^21 from
 * Physical Review B 89, 184102 (2014)
 */
class PolarPFMPhi : public DerivativeParsedMaterialHelper, public ExpressionBuilder
{
public:
  static InputParameters validParams();

  PolarPFMPhi(const InputParameters & parameters);

protected:
  /// coupled interfacial melt order parameter
  EBTerm _upsilon;

  ///@{ interpolation coefficients
  const Real _a_phi;
  const Real _a0;
  ///@}

  /// solid 2 / solid 1 gradient energy coefficient
  const Real _beta21;
};
