/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"

class FFTBufferAux;
template <typename T>
class FFTBufferBase;

template <>
InputParameters validParams<FFTBufferAux>();

/**
 *
 */
class FFTBufferAux : public AuxKernel
{
public:
  FFTBufferAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const FFTBufferBase<Real> & _buffer;
};
