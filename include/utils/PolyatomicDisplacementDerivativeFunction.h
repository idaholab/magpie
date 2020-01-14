/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include "PolyatomicDisplacementFunctionBase.h"
#include "PolyatomicDisplacementFunction.h"

/**
 * Implements the computation of the derivative of the polyatomic displacement functions
 * w.r.t. to composition using Parkin & Coulter's method, compare JNM 101, 1981.
 */
class PolyatomicDisplacementDerivativeFunction : public PolyatomicDisplacementFunctionBase
{
public:
  /// default constructor
  PolyatomicDisplacementDerivativeFunction(std::vector<MyTRIM_NS::Element> polyatomic_material,
                                           nrt_type damage_function_type,
                                           const PolyatomicDisplacementFunction * net_disp,
                                           std::vector<std::vector<Real>> Ecap = {{}});

  virtual std::vector<Real> getRHS(Real energy) override;

  /// computes term 1 in Parkin-Coulter expression nu_k(T - Eb)
  Real
  integralTypeI(Real energy, unsigned int i, unsigned int j, unsigned int l, unsigned int k) const;

  /// computes terms 2 & 3  in Parkin-Coulter expression (nu_i(E-T) & nu_i(E) terms)
  Real
  integralTypeII(Real energy, unsigned int i, unsigned int j, unsigned int l, unsigned int k) const;

  /// computes contribution to derivative computation that do not depend on the variable itself
  Real source(Real energy, unsigned int i, unsigned int j, unsigned int l);

protected:
  /// the source term in the NRT equatons for the derivative requires the solution of the equations themselves
  const PolyatomicDisplacementFunction * _net_displacement_function;
};

#endif // GSL_ENABLED
