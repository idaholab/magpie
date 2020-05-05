/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#pragma once

#include "FFTBufferBase.h"
#include "PerfGraphInterface.h"

#include "fftw3.h"

template <typename T>
class FFTWBufferBase;

/**
 * FFTW specific interleaved data buffer base class
 */
template <typename T>
class FFTWBufferBase : public FFTBufferBase<T>
{
public:
  FFTWBufferBase(const InputParameters & parameters);
  ~FFTWBufferBase();

  // transforms
  void forwardRaw() override;
  void backwardRaw() override;

  // scaling
  Real forwardScale() override { return _scaling; }

protected:
  ///@{ FFTW plans
  fftw_plan _forward_plan;
  fftw_plan _backward_plan;
  ///@}

  ///@{ timers
  PerfID _perf_plan;
  PerfID _perf_fft;
  ///@}

  /// scale factor
  Real _scaling;

  usingFFTBufferBaseMembers;
};

#endif
