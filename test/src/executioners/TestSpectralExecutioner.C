/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "TestSpectralExecutioner.h"
#include "InputParameters.h"

// testing
registerMooseObject("MagpieTestApp", TestSpectralExecutioner);

InputParameters
TestSpectralExecutioner::validParams()
{
  InputParameters params = SpectralExecutionerBase::validParams();
  params.addClassDescription("Spectral executioner for FFT spatial derivatives.");
  return params;
}

TestSpectralExecutioner::TestSpectralExecutioner(const InputParameters & parameters)
  : SpectralExecutionerBase(parameters)
{
}

void
TestSpectralExecutioner::execute()
{
  _time_step = 0;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.advanceState();

  // back and forth test
  auto & c_buffer = getFFTBuffer<Real>("c");
  auto c = c_buffer.realSpace();
  // auto c_tilde = c_buffer.reciprocalSpace();
  c_buffer.forward();

  auto & R_buffer = getFFTBuffer<RealVectorValue>("R");
  auto R = R_buffer.realSpace();
  R_buffer.forward();

  // gradient tests
  auto & u_buffer = getFFTBuffer<Real>("u");
  auto u = u_buffer.realSpace();
  u_buffer.forward();
  auto & v_buffer = getFFTBuffer<Real>("v");
  auto v = v_buffer.realSpace();
  v_buffer.forward();

  auto & grad_u_buffer = getFFTBuffer<RealVectorValue>("grad_u");
  kVectorMultiply(u_buffer, grad_u_buffer);
  auto & grad_v_buffer = getFFTBuffer<RealVectorValue>("grad_v");
  kVectorMultiply(v_buffer, grad_v_buffer);

  _time_step = 1;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
  _fe_problem.advanceState();

  // back and forth test
  c_buffer.backward();
  R_buffer.backward();

  // gradient test
  u_buffer.backward();
  grad_u_buffer.backward();
  v_buffer.backward();
  grad_v_buffer.backward();

  _time_step = 2;
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
}
