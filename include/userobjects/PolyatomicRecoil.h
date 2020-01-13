/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include "ParkinCoulterBase.h"

class PolyatomicRecoil;

class PolyatomicRecoil : public ParkinCoulterBase
{
public:
  static InputParameters validParams();

  PolyatomicRecoil(const InputParameters & parameters);
  void finalize() override;
  void execute() override;

protected:
  virtual void initDamageFunctions() override;
  virtual std::vector<unsigned int> atomicNumbers() const override;
  virtual std::vector<Real> massNumbers() const override;
  virtual std::vector<Real> numberFractions() const override;
  virtual Real maxEnergy() const override { return getParam<Real>("Emax"); }
};

#endif // GSL_ENABLED
