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
  mooseError("SpectralExecutionerBase is an abstract base class");
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
  const Complex I(0.0, 1.0);
  switch (in_buffer.dim())
  {
    case 1:
    {
      const auto & ivec = in_buffer.kTable(0);
      const int ni = grid[0];
      for (int i = 0; i * 2 <= ni; ++i)
        out[i](0) = in[i] * ivec[i] * I;
      return;
    }

    case 2:
    {
      std::size_t index = 0;
      const auto & ivec = in_buffer.kTable(0);
      const auto & jvec = in_buffer.kTable(1);
      const int ni = grid[0];
      const int nj = (grid[1] >> 1) + 1;
      for (int i = 0; i < ni; ++i)
        for (int j = 0; j * 2 < nj; ++j)
        {
          out[index](0) = in[index] * ivec[i] * I;
          out[index](1) = in[index] * jvec[j] * I;
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
            out[index](0) = in[index] * ivec[i] * I;
            out[index](1) = in[index] * jvec[j] * I;
            out[index](2) = in[index] * kvec[k] * I;
            index++;
          }
      return;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}
