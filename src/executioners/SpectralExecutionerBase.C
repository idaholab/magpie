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
  auto & c = getFFTBuffer<Real>("c");
  c.forward();
  auto & R = getFFTBuffer<RealVectorValue>("R");
  R.forward();

  // gradient test
  auto & u = getFFTBuffer<Real>("u");
  u.forward();
  auto & grad_u = getFFTBuffer<RealVectorValue>("grad_u");
  kVectorMultiply(u, grad_u);

  _time_step = 1;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
  _fe_problem.advanceState();

  // back and forth test
  c.backward();
  R.backward();
  R /= 10000.0;
  c /= 10000.0;

  // gradient test
  u.backward();
  grad_u.backward();
  u /= 10000.0;
  grad_u /= 100.0;

  _time_step = 2;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
}

void
SpectralExecutionerBase::kVectorMultiply(const FFTBufferBase<Real> & in,
                                         FFTBufferBase<RealVectorValue> & out) const
{
  mooseAssert(in.size() == out.size(), "Buffer sizes must be equal");
  mooseAssert(in.dim() == out.dim(), "Buffer dimensions must be equal");

  const auto & grid = in.grid();
  switch (in.dim())
  {
    case 1:
    {
      const int ni = grid[0];
      for (int i = 0; i < ni; ++i)
      {
        out[i](0) = in[i] * i;
        out[i](1) = 0.0;
        out[i](2) = 0.0;
      }
      return;
    }

    case 2:
    {
      std::size_t index = 0;
      const int ni = grid[0];
      const int nj = grid[1];
      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
        {
          out[index](0) = in[index] * i;
          out[index](1) = in[index] * j;
          out[index](2) = 0.0;
          index++;
        }
      return;
    }

    case 3:
    {
      std::size_t index = 0;
      const int ni = grid[0];
      const int nj = grid[1];
      const int nk = grid[2];
      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
          for (int k = 0; k < nk; ++k)
          {
            out[index](0) = in[index] * i;
            out[index](1) = in[index] * j;
            out[index](2) = in[index] * k;
            index++;
          }
      return;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}
