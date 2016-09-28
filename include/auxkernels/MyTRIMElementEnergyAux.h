#ifndef MYTRIMELEMENTENERGYAUX_H
#define MYTRIMELEMENTENERGYAUX_H

#include "MyTRIMElementEnergyAccess.h"
#include "AuxKernel.h"

// forward declarations
class MyTRIMElementRun;
class MyTRIMElementEnergyAux;

template<>
InputParameters validParams<MyTRIMElementEnergyAux>();

class MyTRIMElementEnergyAux : public MyTRIMElementEnergyAccess<AuxKernel>
{
public:
  MyTRIMElementEnergyAux(const InputParameters & params);
  virtual ~MyTRIMElementEnergyAux() {}

  virtual Real computeValue();
};

#endif //MYTRIMELEMENTENERGYAUX_H
