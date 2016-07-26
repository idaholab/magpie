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

#ifndef DISCRETEFISSIONPDFTEST_H
#define DISCRETEFISSIONPDFTEST_H

//CPPUnit includes
#include "cppunit/extensions/HelperMacros.h"
#include "DiscreteFissionPKAPDF.h"
#include "MooseRandom.h"

class DiscreteFissionPDFTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DiscreteFissionPDFTest );

  CPPUNIT_TEST( sampleFissionPKA );

  CPPUNIT_TEST_SUITE_END();

public:
  void sampleFissionPKA();
  void setRandomDirection(MyTRIM_NS::IonBase & ion);
};

#endif  // DiscreteFissionPDFTest_H
