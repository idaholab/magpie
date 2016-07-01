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

#include "ReferenceWrapperTest.h"

//Magpie includes
#include <MagpieUtils.h>

CPPUNIT_TEST_SUITE_REGISTRATION( ReferenceWrapperTest );

void
ReferenceWrapperTest::testRefWrapper()
{
  unsigned int l = 0;
  MagpieUtils::reference_wrapper<unsigned int> ref_wrapper_l(l);

  // the simplest of tests
  CPPUNIT_ASSERT( ref_wrapper_l.get() == l );
  //CPPUNIT_ASSERT( ref_wrapper_l == l );

  // make sure that when l changes we pick that change up
  l++;
  CPPUNIT_ASSERT( ref_wrapper_l.get() == l );

  // test the assignment operator
  ref_wrapper_l = 5;
  CPPUNIT_ASSERT( l == 5 );

  // check the += / -= / *= / /= / %=
  ref_wrapper_l += 2;
  CPPUNIT_ASSERT( l == 7 );
  ref_wrapper_l -= 1;
  CPPUNIT_ASSERT( l == 6 );
  ref_wrapper_l *= 3;
  CPPUNIT_ASSERT( l == 18 );
  ref_wrapper_l /= 2;
  CPPUNIT_ASSERT( l == 9 );
  ref_wrapper_l = 13;
  ref_wrapper_l %= 10;
  CPPUNIT_ASSERT( l == 3 );

  // bitwise operators <<=, >>=
  unsigned int l1 = 12;
  ref_wrapper_l = 12;
  l1 <<= 1;
  ref_wrapper_l <<= 1;
  CPPUNIT_ASSERT( l == l1 );
  l1 >>= 2;
  ref_wrapper_l >>= 2;
  CPPUNIT_ASSERT( l == l1 );

  // reassign to a different value
  unsigned int l2 = 123;
  ref_wrapper_l.set(l2);
  ref_wrapper_l += 1;
  CPPUNIT_ASSERT( l2 == 124 );

  // bitwise and &=: 11 & 13 = 9
  unsigned int l3 = 11;
  ref_wrapper_l.set(l3);
  ref_wrapper_l &= 13;
  CPPUNIT_ASSERT( l3 == 9 );

  // bitwise or |=: 9 & 12 = 13
  ref_wrapper_l = 9;
  ref_wrapper_l.set(l3);
  ref_wrapper_l |= 12;
  CPPUNIT_ASSERT( l3 == 13 );

  // bitwise xor ^=: 11 & 13 = 6
  ref_wrapper_l = 11;
  ref_wrapper_l.set(l3);
  ref_wrapper_l ^= 13;
  CPPUNIT_ASSERT( l3 == 6 );
}
