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
 * Gradient energy term in the polar phase field model from Physical Review B 89, 184102 (2014)
 */
class PolarPFMGradient : public DerivativeMaterialInterface<KernelValue>
{
public:
  static InputParameters validParams();

  PolarPFMGradient(const InputParameters & parameters);

protected:
  virtual Real precomputeQpResidual() override;
  virtual Real precomputeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  ///@{ coupled (other) order parameter
  const VariableGradient & _grad_v;
  std::string _v_name;
  unsigned int _v_var;
  ///@}

  ///@{ property derivivatives
  const MaterialProperty<Real> & _dpropdu;
  const MaterialProperty<Real> & _d2propdu2;
  const MaterialProperty<Real> & _d2propdudv;
  ///@}
};
