/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"
#include "MyTRIMRasterizer.h"
#include "mytrim/simconf.h"

class MyTRIMDensityAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MyTRIMDensityAux(const InputParameters & params);

  virtual Real computeValue();

protected:
  const MyTRIMRasterizer & _rasterizer;
  const MyTRIMRasterizer::TrimParameters & _trim_parameters;

  /// number of elements used in the problem
  unsigned int _nvars;

private:
  /// internal TRIM simulation status object
  MyTRIM_NS::SimconfType _simconf;

  /// calculate values only for qp 0 and cache them here
  Real _value_cache;
};
