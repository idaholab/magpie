/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMElementResultAccess.h"
#include "MyTRIMRasterizer.h"
#include "Kernel.h"

class MyTRIMElementSource : public MyTRIMElementResultAccess<Kernel>
{
public:
  static InputParameters validParams();

  MyTRIMElementSource(const InputParameters & params);

protected:
  virtual Real computeQpResidual();

  /// factor for conversion between defect number densities and concentrations
  const Real _prefactor;

  /// Simulation parameters
  const MyTRIMRasterizer::TrimParameters & _trim_parameters;
};
