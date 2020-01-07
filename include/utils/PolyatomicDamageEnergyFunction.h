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

/**
 * Implements the computation of polyatomic displacement functions using Parkin & Coulter's
 * method, compare JNM 101, 1981.
 */
class PolyatomicDamageEnergyFunction : public PolyatomicDisplacementFunctionBase
{
public:
  /// default constructor
  PolyatomicDamageEnergyFunction(std::vector<MyTRIM_NS::Element> polyatomic_material,
                                 nrt_type damage_function_type,
                                 std::vector<std::vector<Real>> Ecap = {{}});

  static int odeRHS(Real energy, const Real disp[], Real f[], void * params);

  /// a getter needed for accessing this pointer in odeRHS
  Real taylorSeriesThreshold() const { return _taylor_series_threshold; }

  Real integralTypeI(Real energy, unsigned int i, unsigned int j);

  Real integralTypeII(Real energy, unsigned int i, unsigned int j);

  /**
   * energy threshold below which the contribution from nu_i(E-T) - nu_i(E) is approximated
   * by nu_i(E-T) - nu_i(E) \approx -T d(nu_i) / dE.
   */
  Real _taylor_series_threshold = 1e2;

protected:
  unsigned int mapIndex(unsigned int i, unsigned int /*j*/, unsigned int /*l*/) const override
  {
    return i;
  };
};

#endif // GSL_ENABLED
