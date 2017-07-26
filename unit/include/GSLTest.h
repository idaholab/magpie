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

#ifndef GSLTEST_H
#define GSLTEST_H

//CPPUnit includes
#include "cppunit/extensions/HelperMacros.h"

class GSLTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( GSLTest );

  CPPUNIT_TEST( integrationTest );

  CPPUNIT_TEST_SUITE_END();

public:
  void integrationTest();

protected:
  static double function(double x, void * params);
};

#endif  // GSLTEST_H
