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

#include "PolyatomicDisplacementFunctionBase.h"

// mytrim includes
#include <mytrim/element.h>

class PolyatomicDisplacementFunction;
class PolyatomicDamageEnergyFunction;
class PolyatomicDisplacementDerivativeFunction;

class ParkinCoulterInterface
{
public:
  ParkinCoulterInterface(const MooseObject * moose_object);
  static InputParameters validParams();

protected:
  /// recomputes the Polyatomic damage functions
  void computeDamageFunctions();

  /// this function is called in computeDamageFunctions and is intended to allow
  /// allocation of _padf &| _padf_derivative
  virtual void initDamageFunctions() = 0;

  ///@{ these functions provide Z, A, and number fraction and must be overidden
  virtual std::vector<unsigned int> atomicNumbers() const = 0;
  virtual std::vector<Real> massNumbers() const = 0;
  virtual std::vector<Real> numberFractions() const = 0;
  ///@}

  /// this function provides the maximum energy to integrate PC eqs to
  virtual Real maxEnergy() const = 0;

  /// computes the polymat object that is used to create PolyatomicDisplacementFunctions
  std::vector<MyTRIM_NS::Element> polyMat() const;

  /// the MooseObject that inherits from this interface
  const MooseObject * _moose_obj;

  /// convenience reference to InputParameters from MooseObject
  const InputParameters & _pars;

  std::vector<std::vector<Real>> _Ecap;

  std::unique_ptr<PolyatomicDisplacementFunctionBase> _padf;
  std::unique_ptr<PolyatomicDisplacementDerivativeFunction> _padf_derivative;
};

#endif // GSL_ENABLED
