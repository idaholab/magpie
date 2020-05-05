/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "SpectralExecutionerBase.h"

// System includes
#include <string>

// Forward declarations
class InputParameters;

/**
 * Executioner for diffusion spectral solver.
 */
class SpectralExecutionerDiffusion : public SpectralExecutionerBase
{
public:
  static InputParameters validParams();

  SpectralExecutionerDiffusion(const InputParameters & parameters);
  /**
   * Algorithm for incremental solution using forward/backward transforms of Green's function.
   */
  virtual void execute() final;
  /**
   * Algorithm for explicit diffusion solve through finite differences.
   */
  virtual void executeExplicit();
  /**
   * Algorithm for non-incremental numerical integration.
   */
  virtual void executeTotalDiffusion();

protected:
  /// Generic diffusion coefficient
  const Real _diff_coeff;

  /// Time step
  const Real _dt;

  /// Number of steps
  const unsigned int _nsteps;

  /// Current time
  Real _t_current;

private:
  /**
   * Helper function to advance diffusion step in an explicit finite difference sense.
   */
  void advanceDiffusionStep(FFTBufferBase<Real> & inOut, const Real diff_coeff);
  /**
   * Helper function to advance diffusion step in a total fashion --no increments.
   */
  void advanceDiffusionStepTotal(const FFTBufferBase<Real> & in,
                                 FFTBufferBase<Real> & out,
                                 const Real D,
                                 const Real time);
  /**
   * Helper function to get the diffusion equation Green's function corresponding to one time step.
   */
  void getGreensFunction(FFTBufferBase<Real> & greens, Real time, const Real D);
  /**
   * Helper function to multiply a buffer by a real number.
   */
  void scalarMultiplyBuffer(FFTBufferBase<Real> & u, const FFTBufferBase<Real> & greens);
};
