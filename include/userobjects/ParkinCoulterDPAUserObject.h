/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include "DPAUserObjectBase.h"
#include "ParkinCoulterInterface.h"

class ParkinCoulterDPAUserObject;

template <>
InputParameters validParams<ParkinCoulterDPAUserObject>();

class ParkinCoulterDPAUserObject : public DPAUserObjectBase, public ParkinCoulterInterface
{
public:
  ParkinCoulterDPAUserObject(const InputParameters & parameters);
  void finalize() override;
  void execute() override;
  void initialSetup() override;

protected:
  virtual void initDamageFunctions() override;
  virtual std::vector<unsigned int> atomicNumbers() const override;
  virtual std::vector<Real> massNumbers() const override;
  virtual std::vector<Real> numberFractions() const override;
  virtual Real maxEnergy() const override { return getMaxEnergy(); }
  virtual Real integralDamageFunction(Real T, unsigned int i, unsigned int j) const override;
  virtual void onCompositionChanged() override;
};

#endif // GSL_ENABLED
