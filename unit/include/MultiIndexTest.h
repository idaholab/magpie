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

#ifndef MULTIINDEXTEST_H
#define MULTIINDEXTEST_H

//CPPUnit includes
#include "cppunit/extensions/HelperMacros.h"

template <class T> class MultiIndex;

class MultiIndexTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MultiIndexTest );

  CPPUNIT_TEST( setUp );
  CPPUNIT_TEST( testIterators );
  CPPUNIT_TEST( dataStoreLoad );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void testIterators();
  void dataStoreLoad();
};

#endif  // MULTIINDEXTEST_H
