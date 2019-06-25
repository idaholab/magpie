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

  static int odeRHS(Real energy, const Real disp[], Real f[], void * params);

  ///@{ some getters needed for accessing this pointer in odeRHS
  Real linearInterpolation(unsigned int i, unsigned int j, unsigned int l, Real energy) const;
  Real linearInterpolation(
      unsigned int i, unsigned int j, unsigned int l, Real energy, unsigned int index) const;
  ///@}

  /// computes term 1 in Parkin-Coulter expression nu_k(T - Eb)
  Real
  integralTypeI(Real energy, unsigned int i, unsigned int j, unsigned int l, unsigned int k) const;

  /// computes terms 2 & 3  in Parkin-Coulter expression (nu_i(E-T) & nu_i(E) terms)
  Real
  integralTypeII(Real energy, unsigned int i, unsigned int j, unsigned int l, unsigned int k) const;

  /// computes contribution to derivative computation that do not depend on the variable itself
  Real source(Real energy, unsigned int i, unsigned int j, unsigned int l);

protected:
  /// maps triple index theta_ijl to single theta_n; l runs faster than j runs faster than i
  unsigned int mapIndex(unsigned int i, unsigned int j, unsigned int l) const
  {
    return i + j * _n_species + l * _n_species * _n_species;
  };

  /// the source term in the NRT equatons for the derivative requires the solution of the equations themselves
  const PolyatomicDisplacementFunction * _net_displacement_function;
};

#endif // PolyatomicDisplacementDerivativeFunction
#endif
