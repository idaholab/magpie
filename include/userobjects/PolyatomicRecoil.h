/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#ifndef POLYATOMICRECOIL_H
#define POLYATOMICRECOIL_H

#include "GeneralUserObject.h"

class PolyatomicRecoil;
class PolyatomicDisplacementFunction;
class PolyatomicDamageEnergyFunction;
class PolyatomicDisplacementDerivativeFunction;

template <>
InputParameters validParams<PolyatomicRecoil>();

class PolyatomicRecoil : public GeneralUserObject
{
public:
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

#endif // PolyatomicRecoil_H
#endif
