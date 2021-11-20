/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"

template <typename T>
class FFTBufferBase;

/**
 *
 */
class FFTBufferAux : public AuxKernel
{
public:
  FFTBufferAux(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeValue() override;

  /// buffer reference as generic parent class
  const UserObject & _buffer;

  ///@{ specific buffer type, only one of those will be non-null
  const FFTBufferBase<Real> * _real_buffer;
  const FFTBufferBase<RealVectorValue> * _realvectorvalue_buffer;
  const FFTBufferBase<RankTwoTensor> * _ranktwotensor_buffer;
  const FFTBufferBase<RankThreeTensor> * _rankthreetensor_buffer;
  const FFTBufferBase<RankFourTensor> * _rankfourtensor_buffer;
  ///@}

  /// component to access
  const std::vector<unsigned int> _component;
};
