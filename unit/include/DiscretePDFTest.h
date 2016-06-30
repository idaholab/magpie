#ifndef DISCRETEPDFTEST_H
#define DISCRETEPDFTEST_H

//CPPUnit includes
#include "cppunit/extensions/HelperMacros.h"
#include "DiscretePKAPDF.h"

class DiscretePDFTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DiscretePDFTest );

  CPPUNIT_TEST( samplePKA );

  CPPUNIT_TEST_SUITE_END();

public:
  void samplePKA();
};

#endif  // DiscretePDFTest_H
