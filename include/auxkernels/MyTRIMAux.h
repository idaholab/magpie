#ifndef MYTRIMAUX_H
#define MYTRIMAUX_H

#include "AuxKernel.h"

// forward declarations
class MyTRIMRun;
class MyTRIMAux;

template<>
InputParameters validParams<MyTRIMAux>();

class MyTRIMAux : public AuxKernel
{
public:
  MyTRIMAux(const InputParameters & params);
  virtual ~MyTRIMAux() {}

  virtual Real computeValue();

protected:
  const MyTRIMRun & _mytrim;
  const unsigned int _ivar;
  const unsigned int _defect;

private:
  Real _value_cache;
};

#endif //MYTRIMAUX_H
