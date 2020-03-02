/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "DerivativeMaterialInterface.h"
#include "KernelValue.h"

/**
 * Bulk energy derivative term in the polar phase field model from Physical Review B 89, 184102
 * (2014)
 */
class PolarPFMDerivative : public DerivativeMaterialInterface<KernelValue>
{
public:
  static InputParameters validParams();

  PolarPFMDerivative(const InputParameters & parameters);

protected:
  virtual Real precomputeQpResidual() override;
  virtual Real precomputeQpJacobian() override;

  const MaterialProperty<Real> & _dpropdu;
  const MaterialProperty<Real> & _d2propdu2;
};
