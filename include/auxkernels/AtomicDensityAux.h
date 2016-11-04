#ifndef ATOMICDENSITYAUX_H
#define ATOMICDENSITYAUX_H

#include "AuxKernel.h"

class AtomicDensityAux;
class MyTRIMRasterizer;

template<>
InputParameters validParams<AtomicDensityAux>();

/**
 * Compute the atomic density at an element.
 */
class AtomicDensityAux : public AuxKernel
{
public:
  AtomicDensityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// Rasterizer object to provide the material data
  const MyTRIMRasterizer & _rasterizer;
};

#endif // ATOMICDENSITYAUX_H
