/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMELEMENTHEATSOURCE_H
#define MYTRIMELEMENTHEATSOURCE_H

#include "MyTRIMElementEnergyAccess.h"
#include "Kernel.h"

// forward declarations
class MyTRIMElementRun;
class MyTRIMElementHeatSource;

template<>
InputParameters validParams<MyTRIMElementHeatSource>();

class MyTRIMElementHeatSource : public MyTRIMElementEnergyAccess<Kernel>
{
public:
  MyTRIMElementHeatSource(const InputParameters & params);

protected:
  virtual Real computeQpResidual();

  // current timestep size
  const Real & _dt;
};

#endif //MYTRIMELEMENTHEATSOURCE_H
