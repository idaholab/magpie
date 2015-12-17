#ifndef SPPARKSAUX_H
#define SPPARKSAUX_H

#include "AuxKernel.h"

// forward declarations
class SPPARKSUserObject;
class SPPARKSAux;

template<>
InputParameters validParams<SPPARKSAux>();

class SPPARKSAux : public AuxKernel
{
public:
  SPPARKSAux(const InputParameters & params);
  virtual ~SPPARKSAux() {}

  virtual Real computeValue();

protected:
  const SPPARKSUserObject & _spparks;
  const unsigned int _var;
  const unsigned int _array;
};

#endif //SPPARKSAUX_H
