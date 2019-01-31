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

// Magpie includes
#include "overlap.hpp"
#include <gtest/gtest.h>

TEST(OverlapTest, tet)
{
  // construct tet
  OVERLAP::vector_t v0{0, 0, 0};
  OVERLAP::vector_t v1{1.2, 0, 0};
  OVERLAP::vector_t v2{0, 1.2, 0};
  OVERLAP::vector_t v3{0, 0, 1.2};
  OVERLAP::Tetrahedron tet = OVERLAP::Tetrahedron{v0, v1, v2, v3};

  // construct sphere
  OVERLAP::Sphere sph(OVERLAP::vector_t{0.0, 0.0, 0.0}, 1.0);

  // answer of 0.283529 is confirmed by Monte-Carlo
  EXPECT_NEAR(OVERLAP::overlap(sph, tet), 0.283529, 1e-5) << "Tet overlap is wrong";
}
