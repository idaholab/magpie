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

  mooseInfo("SpectralExecutionerBase::execute()");

  auto & c = getFFTBuffer<Real>("c");
  c.forward();

  _time_step = 1;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
}
