#ifndef PKAAUX_H
#define PKAAUX_H

#include "AuxKernel.h"

// forward declarations
class MyTRIMRasterizer;
class PKAAux;

template<>
InputParameters validParams<PKAAux>();

class PKAAux : public AuxKernel
{
public:
  PKAAux(const InputParameters & params);
  virtual ~PKAAux() {}

  virtual Real computeValue();

protected:
  const MyTRIMRasterizer & _rasterizer;

};

#endif //PKAAUX_H
