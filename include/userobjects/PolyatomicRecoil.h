/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include "GeneralUserObject.h"

class PolyatomicDisplacementFunctionBase;
class PolyatomicDisplacementDerivativeFunction;

class PolyatomicRecoil : public GeneralUserObject
{
public:
  static InputParameters validParams();

  PolyatomicRecoil(const InputParameters & parameters);

  void execute() override;
  void initialize() override {}
  void finalize() override;

protected:
  std::vector<unsigned int> _atomic_numbers;
  std::vector<Real> _mass_numbers;

  std::unique_ptr<PolyatomicDisplacementFunctionBase> _padf;
  std::unique_ptr<PolyatomicDisplacementDerivativeFunction> _padf_derivative;
};

#endif // GSL_ENABLED
