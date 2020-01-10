/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "PolyatomicDisplacementFunction.h"
#include <gtest/gtest.h>
#include <mytrim/element.h>

TEST(DisplacementFunctionTest, integralDispFunction)
{
  std::vector<unsigned int> Z = {6, 14};
  std::vector<Real> A = {12.0, 28.0};
  std::vector<Real> N = {0.5, 0.5};
  std::vector<Real> threshold = {16.3, 92.6};
  std::vector<MyTRIM_NS::Element> poly_mat;
  for (unsigned int j = 0; j < Z.size(); ++j)
  {
    MyTRIM_NS::Element element;
    element._Z = Z[j];
    element._m = A[j];
    element._t = N[j];
    element._Edisp = threshold[j];
    element._Elbind = 0.0;
    poly_mat.push_back(element);
  }

  PolyatomicDisplacementFunction padf(poly_mat, NET);

  Real energy = padf.minEnergy();
  for (unsigned int j = 0; j < 51; ++j)
  {
    energy += 1;
    padf.advanceDisplacements(energy);
  }
  padf.computeDisplacementFunctionIntegral();
  unsigned int final = padf.nEnergySteps() - 1;

  Real integral = 0;
  for (unsigned int j = 0; j <= final; ++j)
  {
    Real mp = (j == 0 || j == final) ? 0.5 : 1;
    Real energy = padf.energyPoint(j);
    integral += mp * padf.linearInterpolation(energy, 0, 0, 0);
  }
  EXPECT_NEAR(integral,
              padf.linearInterpolationIntegralDamageFunction(padf.energyPoint(final), 0, 0, 0),
              1e-4)
      << "Displacement function integration result is wrong";
}
