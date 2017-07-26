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

#include "GSLTest.h"
#include <cmath>

//GSL includes
#include <gsl/gsl_integration.h>

CPPUNIT_TEST_SUITE_REGISTRATION( GSLTest );

void
GSLTest::integrationTest()
{
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  double result, error;
  double alpha = 1.0;

  gsl_function F;
  F.function = &GSLTest::function;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(result, -4.000000000000085265, 1e-18);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(error, 0.000000000000135447, 1e-18);
  CPPUNIT_ASSERT(w->size == 8);

  gsl_integration_workspace_free (w);
}

// function to integrate
double
GSLTest::function(double x, void * params)
{
  double alpha = *(double *) params;
  double f = std::log(alpha*x) / std::sqrt(x);
  return f;
}
