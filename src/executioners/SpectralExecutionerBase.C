/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "SpectralExecutionerBase.h"
#include "InputParameters.h"

// testing
registerMooseObject("MagpieApp", SpectralExecutionerBase);

InputParameters
SpectralExecutionerBase::validParams()
{
  InputParameters params = Executioner::validParams();
  params.addClassDescription("Executioner for FFT simulations.");
  params.addParam<Real>("time", 0.0, "System time");
  return params;
}

SpectralExecutionerBase::SpectralExecutionerBase(const InputParameters & parameters)
  : Executioner(parameters),
    _system_time(getParam<Real>("time")),
    _time_step(_fe_problem.timeStep()),
    _time(_fe_problem.time()),
    _final_timer(registerTimedSection("final", 1)),
    _fft_problem(dynamic_cast<FFTProblem *>(&_fe_problem))
{

  if (!_fft_problem)
    mooseError("Use Problem/type=FFTProblem with a spectral executioner");
}

void
SpectralExecutionerBase::init()
{
  if (_app.isRecovering())
  {
    _console << "\nCannot recover FFT solves!\nExiting...\n" << std::endl;
    return;
  }

  // checkIntegrity();
  _fe_problem.execute(EXEC_PRE_MULTIAPP_SETUP);
  _fe_problem.initialSetup();
}

void
SpectralExecutionerBase::execute()
{
  _time_step = 0;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.advanceState();

  // back and forth test
  auto & c_buffer = getFFTBuffer<Real>("c");
  auto c = c_buffer.realSpace();
  auto c_tilde = c_buffer.reciprocalSpace();
  c_buffer.forward();

  auto & R_buffer = getFFTBuffer<RealVectorValue>("R");
  auto R = R_buffer.realSpace();
  R_buffer.forward();

  // gradient test
  auto & u_buffer = getFFTBuffer<Real>("u");
  auto u = u_buffer.realSpace();
  u_buffer.forward();

  auto & grad_u_buffer = getFFTBuffer<RealVectorValue>("grad_u");
  auto grad_u = grad_u_buffer.realSpace();
  kVectorMultiply(u_buffer, grad_u_buffer);

  _time_step = 1;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
  _fe_problem.advanceState();

  // back and forth test
  c_buffer.backward();
  R_buffer.backward();
  R /= 10000.0;
  c /= 10000.0;

  // gradient test
  u_buffer.backward();
  grad_u_buffer.backward();

  u /= 10000.0;
  grad_u /= 100.0;

  _time_step = 2;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
}

void
SpectralExecutionerBase::kVectorMultiply(const FFTBufferBase<Real> & in_buffer,
                                         FFTBufferBase<RealVectorValue> & out_buffer) const
{
  mooseAssert(in_buffer.dim() == out_buffer.dim(), "Buffer dimensions must be equal");

  const FFTData<Complex> & in = in_buffer.reciprocalSpace();
  FFTData<ComplexVectorValue> & out = out_buffer.reciprocalSpace();
  mooseAssert(in.size() == out.size(), "Buffer sizes must be equal");

  const auto & grid = in_buffer.grid();
  switch (in_buffer.dim())
  {
    case 1:
    {
      const auto & ivec = in_buffer.kTable(0);
      const int ni = grid[0];
      for (int i = 0; i * 2 <= ni; ++i)
        out[i](0) = in[i] * ivec[i];
      return;
    }

    case 2:
    {
      std::size_t index = 0;
      const auto & ivec = in_buffer.kTable(0);
      const auto & jvec = in_buffer.kTable(1);
      const int ni = grid[0];
      const int nj = grid[1];
      for (int i = 0; i < ni; ++i)
        for (int j = 0; j * 2 <= nj; ++j)
        {
          out[index](0) = in[index] * ivec[i];
          out[index](1) = in[index] * jvec[j];
          index++;
        }
      return;
    }

    case 3:
    {
      std::size_t index = 0;
      const auto & ivec = in_buffer.kTable(0);
      const auto & jvec = in_buffer.kTable(1);
      const auto & kvec = in_buffer.kTable(2);
      const int ni = grid[0];
      const int nj = grid[1];
      const int nk = grid[2];
      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
          for (int k = 0; k * 2 <= nk; ++k)
          {
            out[index](0) = in[index] * ivec[i];
            out[index](1) = in[index] * jvec[j];
            out[index](2) = in[index] * kvec[k];
            index++;
          }
      return;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}
