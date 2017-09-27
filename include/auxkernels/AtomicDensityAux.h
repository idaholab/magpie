/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

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

  /// volume scale factor to go from nm^3 to the mesh units selected in the rasterizer
  const Real _volume_scale;
};

#endif // ATOMICDENSITYAUX_H
