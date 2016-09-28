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
};

#endif //MYTRIMELEMENTHEATSOURCE_H
