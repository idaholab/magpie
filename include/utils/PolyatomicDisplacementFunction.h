/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#ifndef POLYATOMICDISPLACEMENTFUNCTION_H
#define POLYATOMICDISPLACEMENTFUNCTION_H

#include "PolyatomicDisplacementFunctionBase.h"

/**
 * Implements the computation of polyatomic displacement functions using Parkin & Coulter's
 * method, compare JNM 101, 1981.
 */
class PolyatomicDisplacementFunction : public PolyatomicDisplacementFunctionBase
{
public:
  /// default constructor
  PolyatomicDisplacementFunction(std::vector<MyTRIM_NS::Element> polyatomic_material,
                                 nrt_type damage_function_type,
                                 std::vector<std::vector<Real>> Ecap = {{}});

  static int odeRHS(Real energy, const Real disp[], Real f[], void * params);

  ///@{ some getters needed for accessing this pointer in odeRHS
  Real linearInterpolation(unsigned int i, unsigned int j, Real energy) const;
  Real linearInterpolation(unsigned int i, unsigned int j, Real energy, unsigned int index) const;
  ///@}

  /// computes term 1 in Parkin-Coulter expression nu_k(T - Eb)
  Real integralTypeI(Real energy, unsigned int i, unsigned int j, unsigned int k) const;

  /// computes terms 2 & 3  in Parkin-Coulter expression (nu_i(E-T) & nu_i(E) terms)
  Real integralTypeII(Real energy, unsigned int i, unsigned int j, unsigned int k) const;

protected:
  /// maps double index nu_ij to single index nu_n; j runs faster than i
  unsigned int mapIndex(unsigned int i, unsigned int j) const { return i + j * _n_species; };

  /// is the total damage function computed
  bool _total_displacement_function;
};

#endif // PolyatomicDisplacementFunction
#endif
