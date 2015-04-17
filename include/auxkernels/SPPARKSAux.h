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
  SPPARKSAux(const std::string & name, InputParameters params);
  virtual ~SPPARKSAux() {}

  virtual Real computeValue();

protected:
  const SPPARKSUserObject & _spparks;
  const int _ivar;
  const int _dvar;
};

#endif
