#ifndef MYTRIMRESULTAUX_H
#define MYTRIMRESULTAUX_H

#include "MyTRIMResultAccess.h"
#include "AuxKernel.h"

// forward declarations
class MyTRIMRun;
class MyTRIMResultAux;

template<>
InputParameters validParams<MyTRIMResultAux>();

class MyTRIMResultAux : public MyTRIMResultAccess<AuxKernel>
{
public:
  MyTRIMResultAux(const InputParameters & params);
  virtual ~MyTRIMResultAux() {}

  virtual Real computeValue();
};

#endif //MYTRIMRESULTAUX_H
