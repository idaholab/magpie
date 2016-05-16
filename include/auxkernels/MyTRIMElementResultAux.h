#ifndef MYTRIMRESULTAUX_H
#define MYTRIMRESULTAUX_H

#include "MyTRIMResultAccess.h"
#include "AuxKernel.h"

// forward declarations
class MyTRIMElementRun;
class MyTRIMElementResultAux;

template<>
InputParameters validParams<MyTRIMElementResultAux>();

class MyTRIMElementResultAux : public MyTRIMResultAccess<AuxKernel>
{
public:
  MyTRIMElementResultAux(const InputParameters & params);
  virtual ~MyTRIMElementResultAux() {}

  virtual Real computeValue();
};

#endif //MYTRIMRESULTAUX_H
