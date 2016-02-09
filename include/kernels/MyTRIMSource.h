#ifndef MYTRIMSOURCE_H
#define MYTRIMSOURCE_H

#include "MyTRIMResultAccess.h"
#include "Kernel.h"

// forward declarations
class MyTRIMRun;
class MyTRIMSource;

template<>
InputParameters validParams<MyTRIMSource>();

class MyTRIMSource : public MyTRIMResultAccess<Kernel>
{
public:
  MyTRIMSource(const InputParameters & params);

protected:
  virtual Real computeQpResidual();

  const Real _prefactor;
};

#endif //MYTRIMAUX_H
