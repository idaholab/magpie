#ifndef MYTRIMSOURCE_H
#define MYTRIMSOURCE_H

#include "MyTRIMResultAccess.h"
#include "Kernel.h"

// forward declarations
class MyTRIMElementRun;
class MyTRIMElementSource;

template<>
InputParameters validParams<MyTRIMElementSource>();

class MyTRIMElementSource : public MyTRIMResultAccess<Kernel>
{
public:
  MyTRIMElementSource(const InputParameters & params);

protected:
  virtual Real computeQpResidual();

  const Real _prefactor;
};

#endif //MYTRIMAUX_H
