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

#ifndef REFERENCEWRAPPERTEST_H
#define REFERENCEWRAPPERTEST_H

//CPPUnit includes
#include "cppunit/extensions/HelperMacros.h"

template <class T> class reference_wrapper;

class ReferenceWrapperTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ReferenceWrapperTest );

  CPPUNIT_TEST( testRefWrapper );
  CPPUNIT_TEST_SUITE_END();

public:
  void testRefWrapper();
};

#endif  // REFERENCEWRAPPERTEST_H
