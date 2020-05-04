/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "Executioner.h"
#include "FFTWBufferBase.h"
#include "FFTProblem.h"

// System includes
#include <string>

// Forward declarations
class InputParameters;

/**
 * FFT Executioner base class.
 */
class SpectralExecutionerBase : public Executioner
{
public:
  static InputParameters validParams();

  SpectralExecutionerBase(const InputParameters & parameters);

  virtual void init() override;
  virtual void execute() override;
  virtual bool lastSolveConverged() const override { return true; }

protected:
  /// obtain a non-const reference to an FFT buffer
  template <typename T>
  FFTBufferBase<T> & getFFTBuffer(const std::string & name);

  /// multiply a scalar buffer by its corresponding k-vector field int a vector buffer
  void kVectorMultiply(const FFTBufferBase<Real> & in, FFTBufferBase<RealVectorValue> & out) const;

  Real _system_time;
  int & _time_step;
  Real & _time;

  PerfID _final_timer;

  FFTProblem * _fft_problem;
};

template <typename T>
FFTBufferBase<T> &
SpectralExecutionerBase::getFFTBuffer(const std::string & name)
{
  return _fft_problem->getFFTBuffer<T>(name);
}
