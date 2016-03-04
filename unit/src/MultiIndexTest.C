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

#include "MultiIndexTest.h"

//Magpie includes
#include "MultiIndex.h"

//libMesh include
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

CPPUNIT_TEST_SUITE_REGISTRATION( MultiIndexTest );

void
MultiIndexTest::setUp()
{
  // Construct a 3 indexed objects of Reals
  std::vector<unsigned int> shape;
  shape.resize(3);
  shape[0] = 3;
  shape[1] = 2;
  shape[2] = 4;
  MagpieUtils::MultiIndex<Real> mindex1 = MagpieUtils::MultiIndex<Real>(shape);

  // check dimension
  CPPUNIT_ASSERT( mindex1.dimension() == 3 );

  // check shape
  std::vector<unsigned int> shape2 = mindex1.shape();
  for (unsigned int j = 0; j < 3; ++j)
    CPPUNIT_ASSERT( shape2[j] == shape[j] );

  // check size operator
  for (unsigned int j = 0; j < 3; ++j)
    CPPUNIT_ASSERT( mindex1.size(j) == shape[j] );

  // parenthesis operator
  std::vector<unsigned int> index;
  index.resize(3);
  for (unsigned int j0 = 0; j0 < shape2[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape2[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape2[2]; ++j2)
      {
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        mindex1(index) = j0 + 10.0 * j1 + 100.0 * j2;
      }

  // check the parenthesis operator but this time reverse loop order
  for (unsigned int j2 = 0; j2 < shape2[2]; ++j2)
    for (unsigned int j0 = 0; j0 < shape2[0]; ++j0)
      for (unsigned int j1 = 0; j1 < shape2[1]; ++j1)
      {
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        CPPUNIT_ASSERT( mindex1(index) == j0 + 10.0 * j1 + 100.0 * j2 );
      }

}
