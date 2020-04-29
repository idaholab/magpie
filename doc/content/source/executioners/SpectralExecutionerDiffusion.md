# SpectralExecutionerDiffusion

!syntax description /Executioner/SpectralExecutionerDiffusion

## Diffusion spectral solve

`SpectralExecutionerDiffusion` inherits from `SpectralExecutionerBase` and implements problem-specific methods. For example, it computes the diffusion equation's Green's function for a time step `_dt` through the method `getGreensFunction`. This function is used in the `execute` method to explicitly obtain the diffusion solution step by step. The Green's function is transformed to the reciprocal space where it is point-wise multiplied by the diffused variable in the Fourier. The transformation of its result into the real space carries the proper grid-related scale factor through the `FFTBufferBase` method `backward`.

As shown below, this diffusion execution does not need to iterate due to its linear nature. Spectral executioners are, however, not restricted to linear solves.

```cpp
void
SpectralExecutionerDiffusion::execute()
{
  unsigned int thisStep = 0;

  _time_step = thisStep;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.advanceState();

  auto & u_buffer = getFFTBuffer<Real>("u");
  auto & greens = getFFTBuffer<Real>("greens_buffer");

  // Get Greens function of
  getGreensFunction(greens, _dt, _diff_coeff);
  greens.forwardRaw();

  for (unsigned int step_no = 0; step_no < _nsteps; step_no++)
  {
    u_buffer.forwardRaw();

    u_buffer.reciprocalSpace() *= greens.reciprocalSpace();
    u_buffer.backward();

    // End of diffusion computations
    thisStep++;
    _t_current += _dt;
    _time_step = thisStep;

    _fe_problem.execute(EXEC_FINAL);
    _time = _t_current;
    Moose::out << "_t_current: " << _t_current << ". \n";
    _fe_problem.outputStep(EXEC_FINAL);

    if (step_no != _nsteps - 1)
      _fe_problem.advanceState();
}
}
```

!syntax parameters /Executioner/SpectralExecutionerDiffusion

!syntax inputs /Executioner/SpectralExecutionerDiffusion

!syntax children /Executioner/SpectralExecutionerDiffusion

!bibtex bibliography
