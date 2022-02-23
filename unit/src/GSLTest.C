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
#ifdef GSL_ENABLED

#include <cmath>
#include <gtest/gtest.h>

// GSL includes
#include <gsl/gsl_integration.h>

// function to integrate
double
function(double x, void * params)
{
  double alpha = *(double *)params;
  double f = std::log(alpha * x) / std::sqrt(x);
  return f;
}

TEST(GSLTest, integrationTest)
{
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  double result, error;
  double alpha = 1.0;

  gsl_function F;
  F.function = &function;
  F.params = &alpha;

  gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);

  EXPECT_NEAR(result, -4.0, 1e-12) << "Integration result is wrong";
  EXPECT_NEAR(error, 0.0, 1e-12) << "Integration error is wrong";
  EXPECT_EQ(w->size, 8);

  gsl_integration_workspace_free(w);
}

#endif
