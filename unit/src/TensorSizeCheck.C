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

#include <cmath>
#include <gtest/gtest.h>

#include "libmesh/vector_value.h"

#include "MooseTypes.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

TEST(TensorSizeCheck, vector) { EXPECT_EQ(sizeof(RealVectorValue), LIBMESH_DIM * sizeof(Real)); }

TEST(TensorSizeCheck, rankTwo)
{
  EXPECT_EQ(sizeof(RankTwoTensor), LIBMESH_DIM * LIBMESH_DIM * sizeof(Real));
}

TEST(TensorSizeCheck, rankThree)
{
  EXPECT_EQ(sizeof(RankThreeTensor), LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM * sizeof(Real));
}

TEST(TensorSizeCheck, rankFour)
{
  EXPECT_EQ(sizeof(RankFourTensor),
            LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM * sizeof(Real));
}
