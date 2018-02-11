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
#include "DiscreteFissionPKAPDF.h"
#include "MooseRandom.h"

#include <gtest/gtest.h>

class DiscreteFissionPDFTest : public ::testing::Test
{
protected:
  void setRandomDirection(MyTRIM_NS::IonBase & ion)
  {
    Real nsq, x1, x2;

    // Marsaglia's method for uniformly sampling the surface of the sphere
    do
    {
      x1 = 2 * MooseRandom::rand() - 1.0;
      x2 = 2 * MooseRandom::rand() - 1.0;
      nsq = x1 * x1 + x2 * x2;
    } while (nsq >= 1);

    // construct normalized direction vector
    ion._dir(0) = 2.0 * x1 * std::sqrt(1.0 - nsq);
    ion._dir(1) = 2.0 * x2 * std::sqrt(1.0 - nsq);
    ion._dir(2) = 1.0 - 2.0 * nsq;
  }
};

#endif  // DiscreteFissionPDFTest_H
