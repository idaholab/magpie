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
  MultiIndex<Real> mindex1 = MultiIndex<Real>(shape);

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

  // check the two-argument constructor
  std::vector<Real> data;
  data.resize(mindex1.nEntries());
  unsigned int p = 0;
  for (unsigned int j0 = 0; j0 < shape2[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape2[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape2[2]; ++j2)
      {
        data[p] = j0 - 3.0 * j1 + 100.0 * j2;
        p++;
      }

  MultiIndex<Real> mindex2 = MultiIndex<Real>(shape, data);
  for (unsigned int j0 = 0; j0 < shape2[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape2[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape2[2]; ++j2)
      {
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        CPPUNIT_ASSERT( mindex2(index) == j0 - 3.0 * j1 + 100.0 * j2 );
      }

  // test resize increasing and descreasing in different dimension simultaneously
  std::vector<unsigned int> new_shape;
  new_shape.resize(3);
  new_shape[0] = 2;
  new_shape[1] = 3;
  new_shape[2] = 5;
  mindex2.resize(new_shape);
  for (unsigned int j0 = 0; j0 < new_shape[0]; ++j0)
    for (unsigned int j1 = 0; j1 < new_shape[1]; ++j1)
      for (unsigned int j2 = 0; j2 < new_shape[2]; ++j2)
      {
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        if (j0 < shape[0] && j1 < shape[1] && j2 < shape[2])
          CPPUNIT_ASSERT( mindex2(index) == j0 - 3.0 * j1 + 100.0 * j2 );
        else
          CPPUNIT_ASSERT( mindex2(index) == 0.0);
      }

      // test assign increasing and descreasing in different dimension simultaneously
      index.resize(4);
      new_shape.resize(4);
      new_shape[0] = 4;
      new_shape[1] = 2;
      new_shape[2] = 1;
      new_shape[3] = 8;
      mindex2.assign(new_shape, 1.0);
      for (unsigned int j0 = 0; j0 < new_shape[0]; ++j0)
        for (unsigned int j1 = 0; j1 < new_shape[1]; ++j1)
          for (unsigned int j2 = 0; j2 < new_shape[2]; ++j2)
            for (unsigned int j3 = 0; j3 < new_shape[3]; ++j3)
            {
              index[0] = j0;
              index[1] = j1;
              index[2] = j2;
              index[3] = j3;
              CPPUNIT_ASSERT( mindex2(index) == 1.0);
            }

}

void
MultiIndexTest::testIterators()
{
  // Create empty MultiIndex object
  std::vector<unsigned int> shape;
  shape.resize(3);
  shape[0] = 3;
  shape[1] = 2;
  shape[2] = 4;
  MultiIndex<Real> mindex = MultiIndex<Real>(shape);

  // set the data by using an iterator
  MultiIndex<Real>::iterator it = mindex.begin();
  MultiIndex<Real>::iterator it_end = mindex.end();
  for (; it != it_end ; ++it)
  {
    std::vector<unsigned int> indices = it.indices();
    *it = indices[0] - 3.0 * indices[1] + 100.0 * indices[2];
  }

  // check the values
  for (unsigned int j0 = 0; j0 < shape[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape[2]; ++j2)
      {
        std::vector<unsigned int> index(3);
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        CPPUNIT_ASSERT( mindex(index) == j0 - 3.0 * j1 + 100.0 * j2 );
      }

  // check the indices function
  it = mindex.begin();
  for (unsigned int j0 = 0; j0 < shape[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape[2]; ++j2)
      {
        std::vector<unsigned int> index(3);
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        CPPUNIT_ASSERT( index[0] == it.index(0) && index[1] == it.index(1) && index[2] == it.index(2));
        ++it;
      }
}

void
MultiIndexTest::dataStoreLoad()
{
  // Create empty MultiIndex object
  std::vector<unsigned int> shape;
  shape.resize(3);
  shape[0] = 3;
  shape[1] = 2;
  shape[2] = 4;
  MultiIndex<Real> mindex = MultiIndex<Real>(shape);

  // set the data by using an iterator
  MultiIndex<Real>::iterator it = mindex.begin();
  MultiIndex<Real>::iterator it_end = mindex.end();
  for (; it != it_end ; ++it)
  {
    std::vector<unsigned int> indices = it.indices();
    *it = indices[0] - 3.0 * indices[1] + 100.0 * indices[2];
  }

  // Serialize
  std::ostringstream oss;
  dataStore(oss, mindex, this);

  // Pour oss into iss
  std::string ostring = oss.str();
  std::istringstream iss(ostring);

  // Clear data structure to avoid false positives and then
  // read data
  for (it = mindex.begin(); it != it_end ; ++it)
    *it = 0.0;

  dataLoad(iss, mindex, this);

  // check the values
  for (unsigned int j0 = 0; j0 < shape[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape[2]; ++j2)
      {
        std::vector<unsigned int> index(3);
        index[0] = j0;
        index[1] = j1;
        index[2] = j2;
        CPPUNIT_ASSERT( mindex(index) == j0 - 3.0 * j1 + 100.0 * j2 );
      }
}
