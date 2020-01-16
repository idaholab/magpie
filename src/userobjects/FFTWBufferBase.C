/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#include "FFTWBufferBase.h"
#include "MyTRIMMesh.h"
#include "MooseTypes.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

template <typename T>
FFTWBufferBase<T>::FFTWBufferBase(const InputParameters & parameters)
  : FFTBufferBase<T>(parameters), _forward_plan(nullptr), _backward_plan(nullptr)
{
  // create plans
}

template <typename T>
FFTWBufferBase<T>::~FFTWBufferBase()
{
  // destroy FFTW plans
  fftw_destroy_plan(_forward_plan);
  fftw_destroy_plan(_backward_plan);
}

template <typename T>
void
FFTWBufferBase<T>::forward()
{
  // execute plan
  fftw_execute(_forward_plan);
}

template <typename T>
void
FFTWBufferBase<T>::backward()
{
  // execute plan
  fftw_execute(_backward_plan);
}

// explicit instantiation and registration
#define FFTWBufferInstance(T)                                                                      \
  template class FFTWBufferBase<T>;                                                                \
  using T##FFTWBuffer = FFTWBufferBase<T>;                                                         \
  registerMooseObject("MagpieApp", T##FFTWBuffer)

FFTWBufferInstance(Real);
FFTWBufferInstance(RealVectorValue);
FFTWBufferInstance(RankTwoTensor);
FFTWBufferInstance(RankThreeTensor);
FFTWBufferInstance(RankFourTensor);

#endif
