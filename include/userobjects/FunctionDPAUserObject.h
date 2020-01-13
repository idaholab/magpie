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
#include "LinearInterpolation.h"

class FunctionDPAUserObject;

template <>
InputParameters validParams<FunctionDPAUserObject>();

class FunctionDPAUserObject : public DPAUserObjectBase
{
public:
  FunctionDPAUserObject(const InputParameters & parameters);
  void finalize() override;
  void execute() override;
  void initialSetup() override;
  void initialize() override {}

protected:
  virtual Real integralDamageFunction(Real T, unsigned int i, unsigned int j) const override;
  virtual void onCompositionChanged() override;

  /// the maximum energy step size used for interpolation and integration of integral damage function
  Real _max_delta_E;

  /// the damage functions are provided by MOOSE functions
  std::vector<std::vector<const Function *>> _damage_functions;

  /// stores the integral damage functions computed from input as LinearInterpolation objects
  std::vector<std::vector<LinearInterpolation>> _integral_damage_functions;
};

#endif
