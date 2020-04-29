# SpectralExecutionerBase

!syntax description /Executioner/SpectralExecutionerBase

## Spectral executioners

In general, executioners determine the type of solve for finite element problems. In spectral solvers, the executioner is also responsible for **actually** solving the system. `SpectralExecutionerBase` is a non-pure base spectral executioners which holds a pointer to an `FFTProblem` and initializes an `FEProblem`.  

The executioner has helper methods such as `kVectorMultiply` which is called to obtain  spatial derivatives in the reciprocal which are then backward transformed into the real space.

The spectral solve application developer will need to inherit from this class and provide appropriate explicit or implicit algorithm to obtain each step's solution.

## ::execute for computing spatial derivatives

The `execute` method handles the numerical time stepping, the forward and inverse discrete Fourier transforms through `FFTWBufferBase` and   

```cpp
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
```

!syntax parameters /Executioner/SpectralExecutionerBase

!syntax inputs /Executioner/SpectralExecutionerBase

!syntax children /Executioner/SpectralExecutionerBase

!bibtex bibliography
