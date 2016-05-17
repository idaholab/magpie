#ifndef MYTRIMRESULTAUX_H
#define MYTRIMRESULTAUX_H

#include "MyTRIMElementResultAccess.h"
#include "AuxKernel.h"

// forward declarations
class MyTRIMElementRun;
class MyTRIMElementResultAux;

template<>
InputParameters validParams<MyTRIMElementResultAux>();

class MyTRIMElementResultAux : public MyTRIMElementResultAccess<AuxKernel>
{
public:
  MyTRIMElementResultAux(const InputParameters & params);
  virtual ~MyTRIMElementResultAux() {}

  virtual Real computeValue();
};

#endif //MYTRIMRESULTAUX_H
