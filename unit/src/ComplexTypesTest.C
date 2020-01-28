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
#include <cmath>
#include <gtest/gtest.h>

#include "ComplexTypes.h"

TEST(ComplexTypesTest, rankTwo)
{
  Complex c_1;
  ComplexType<Real>::type c_2;
  auto c = c_1 + c_2;

  ComplexVectorValue v_1;
  ComplexType<RealVectorValue>::type v_2;
  auto v = v_1 + v_2;

  ComplexRankTwoTensor r2_1;
  ComplexType<RankTwoTensor>::type r2_2;
  auto r2 = r2_1 + r2_2;

  ComplexRankThreeTensor r3_1;
  ComplexType<RankThreeTensor>::type r3_2;
  auto r3 = r3_1 + r3_2;
}
