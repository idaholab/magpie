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
  : FFTBufferBase<T>(parameters),
    _perf_plan(this->registerTimedSection("fftw_plan_r2r", 2)),
    _perf_fft(this->registerTimedSection("fftw_execute", 2))
{
  // create plans
  std::vector<fftw_r2r_kind> kind(_dim, FFTW_R2HC);
  {
    TIME_SECTION(_perf_plan);
    _forward_plan = fftw_plan_r2r(_dim, _grid.data(), _start, _start, kind.data(), FFTW_ESTIMATE);
    _backward_plan = fftw_plan_r2r(_dim, _grid.data(), _start, _start, kind.data(), FFTW_ESTIMATE);
  }
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
  mooseInfo("FFTWBufferBase<T>::forward()");

  // execute plan
  {
    TIME_SECTION(_perf_fft);
    fftw_execute(_forward_plan);
  }
}

template <typename T>
void
FFTWBufferBase<T>::backward()
{
  // execute plan
  {
    TIME_SECTION(_perf_fft);
    fftw_execute(_backward_plan);
  }
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
