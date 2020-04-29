/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "SpectralExecutionerDiffusion.h"
#include "InputParameters.h"

registerMooseObject("MagpieApp", SpectralExecutionerDiffusion);

InputParameters
SpectralExecutionerDiffusion::validParams()
{
  InputParameters params = SpectralExecutionerBase::validParams();
  params.addClassDescription("Executioner for spectral solve of diffusion equation");
  params.addParam<Real>("diffusion_coefficient", 1.0, "Diffusion coefficient");
  params.addParam<Real>("time_step", 1.0, "Time step for ODE integration");
  params.addParam<unsigned int>("number_steps", 0.0, "Time step for ODE integration");
  return params;
}

SpectralExecutionerDiffusion::SpectralExecutionerDiffusion(const InputParameters & parameters)
  : SpectralExecutionerBase(parameters),
    _diff_coeff(getParam<Real>("diffusion_coefficient")),
    _dt(getParam<Real>("time_step")),
    _nsteps(getParam<unsigned int>("number_steps"))
{
  _t_current = 0.0;
}

void
SpectralExecutionerDiffusion::executeExplicit()
{
  unsigned int thisStep = 0;

  _time_step = thisStep;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.advanceState();

  auto & u_buffer = getFFTBuffer<Real>("u");
  auto u = u_buffer.realSpace();

  for (unsigned int step_no = 0; step_no < _nsteps; step_no++)
  {
    // Forward transform to get \hat{u}_{k}
    u_buffer.forwardRaw();

    advanceDiffusionStep(u_buffer, _diff_coeff);

    u_buffer.backward();

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

void
SpectralExecutionerDiffusion::getGreensFunction(FFTBufferBase<Real> & greens,
                                                Real time,
                                                const Real D)
{
  Real accGreens = 0.0;

  const Point & box_size = greens.getBoxSize();
  const Point & min_corner = greens.getMinCorner();
  const Point & max_corner = greens.getMaxCorner();

  const auto & grid = greens.grid();
  switch (greens.dim())
  {
    case 1:
    {
      mooseError("Error: Dimension 1 not implemented yet.");
      break;
    }

    case 2:
    {
      std::size_t index = 0;

      const int ni = grid[0];
      const int nj = grid[1]; // Real space.

      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
        {

          Real x1 = Real(i) / Real(ni) * box_size(0);
          Real y1 = Real(j) / Real(nj) * box_size(1);

          Real x2 = box_size(0) - x1;
          Real y2 = box_size(1) - y1;

          Moose::out << "box size: " << box_size(0) << ", " << box_size(1) << "\n";
          Moose::out << "number of cells: " << ni << ", " << nj << "\n";

          if (time == 0)
            mooseError("Greens function undefined at t = 0s in the FFT solver.");
          else
          {
            Real sum = std::exp(-(x1 * x1 + y1 * y1) / 4 / D / time) +
                       std::exp(-(x1 * x1 + y2 * y2) / 4 / D / time) +
                       std::exp(-(x2 * x2 + y1 * y1) / 4 / D / time) +
                       std::exp(-(x2 * x2 + y2 * y2) / 4 / D / time);

            if (sum == 0 && i != 0 && j != 0)
              mooseError("Precision lost in the analytical integration of heat equation");

            Real correction_factor = (Real(ni) / box_size(0)) * (Real(nj) / box_size(1));

            greens.realSpace()[index] =
                (1.0 / (4 * libMesh::pi * D * time)) * sum / correction_factor;
          }
          accGreens += greens.realSpace()[index];
          index++;
        }
      mooseInfo("Accumulated value of Greens function: ", accGreens, "\n");
      break;
    }

    case 3:
    {
      mooseError("Error: Dimension 3 not implemented yet.");
      break;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}
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
    // scalarMultiplyBuffer(u_buffer, greens);
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

void
SpectralExecutionerDiffusion::executeTotalDiffusion()
{
  unsigned int thisStep = 0;

  _time_step = thisStep;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.advanceState();

  FFTBufferBase<Real> & u_buffer = getFFTBuffer<Real>("u");
  FFTBufferBase<Real> & u_init_buffer = getFFTBuffer<Real>("u_initial");

  u_init_buffer.forwardRaw();
  // auto &greens = getFFTBuffer<Real>("greens_buffer");

  for (unsigned int step_no = 0; step_no < _nsteps; step_no++)
  {
    // Forward transform to get \hat{u}_{k}
    u_buffer.forwardRaw();

    advanceDiffusionStepTotal(u_init_buffer, u_buffer, _diff_coeff, _time);

    u_buffer.backward();

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

void
SpectralExecutionerDiffusion::scalarMultiplyBuffer(FFTBufferBase<Real> & u,
                                                   const FFTBufferBase<Real> & greens)
{
  const auto & grid = u.grid();

  switch (u.dim())
  {
    case 1:
    {
      mooseError("Error: Dimension 1 not implemented yet.");
      return;
    }

    case 2:
    {
      std::size_t index = 0;
      const int ni = grid[0];
      const int nj = (grid[1] >> 1) + 1;

      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
        {
          u.reciprocalSpace()[index] = greens.reciprocalSpace()[index] * u.reciprocalSpace()[index];
          index++;
        }
      return;
    }

    case 3:
    {
      mooseError("Error: Dimension 3 not implemented yet.");
      return;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}

void
SpectralExecutionerDiffusion::advanceDiffusionStepTotal(const FFTBufferBase<Real> & u_initial,
                                                        FFTBufferBase<Real> & u,
                                                        const Real D,
                                                        const Real time)
{
  const auto & grid = u.grid();
  switch (u.dim())
  {
    case 1:
    {
      mooseError("Error: Dimension 1 not implemented yet.");
      return;
    }

    case 2:
    {
      std::size_t index = 0;
      const int ni = grid[0];
      const int nj = (grid[1] >> 1) + 1;

      const auto & ivec = u.kTable(0);
      const auto & jvec = u.kTable(1);

      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
        {
          u.reciprocalSpace()[index] =
              u_initial.reciprocalSpace()[index] *
              std::exp(-D * (ivec[i] * ivec[i] + jvec[j] * jvec[j]) * time);
          index++;
        }
      return;
    }

    case 3:
    {
      mooseError("Error: Dimension 3 not implemented yet.");
      return;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}

void
SpectralExecutionerDiffusion::advanceDiffusionStep(FFTBufferBase<Real> & inOut, const Real D)
{
  const Real sizePixel = 1.0; // Assumed same size in all directions

  const auto & grid = inOut.grid();
  switch (inOut.dim())
  {
    case 1:
    {
      mooseError("Error: Dimension 1 not implemented yet.");
      return;
    }

    case 2:
    {
      std::size_t index = 0;
      const int ni = grid[0];
      const int nj = (grid[1] >> 1) + 1;

      for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
        {
          Real r1 = D * _dt / (sizePixel * sizePixel);
          Real r2 = D * _dt / (sizePixel * sizePixel);

          Real hk1 = (1 - 2 * r1 * (1 - std::cos(libMesh::pi * i / ((ni) / 2))));
          Real hk2 = (1 - 2 * r2 * (1 - std::cos(libMesh::pi * j / ((nj) / 2))));
          inOut.reciprocalSpace()[index] *= hk1 * hk2;
          index++;
        }
      return;
    }

    case 3:
    {
      mooseError("Error: Dimension 3 not implemented yet.");
      return;
    }

    default:
      mooseError("Invalid buffer dimension");
  }
}
