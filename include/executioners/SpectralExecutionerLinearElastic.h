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
class SpectralExecutionerLinearElastic : public SpectralExecutionerBase
{
public:
  static InputParameters validParams();

  SpectralExecutionerLinearElastic(const InputParameters & parameters);
  /**
   * Algorithm for incremental solution using forward/backward transforms of Green's function.
   */
  virtual void execute() final;

protected:
  /// Time step
  const Real _dt;

  /// Number of steps
  const unsigned int _nsteps;

  /// First parameter
  const Real _young_modulus;

  /// Second parameter
  const Real _poisson_ratio;

  /// Current time
  Real _t_current;

  /// Initial homogeneous shear deformation
  const Real _initial_shear_strain;

  /// Initial strain tensor
  RankTwoTensor _initial_strain_tensor;

  /// Average factor to obtain homogeneous material
  const Real _average_factor;

private:
  /**
   * Helper function to get the diffusion equation Green's function corresponding to one time step.
   */
  void getGreensFunction(FFTBufferBase<RankFourTensor> & gamma_hat,
                         FFTBufferBase<Real> & ratio_buffer,
                         const RankFourTensor & elasticity_tensor);

  /**
   * Helper function to get the initial stress from strain and tensor of elastic coefficients.
   */
  FFTBufferBase<RankTwoTensor> & getInitialStress(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                                                  FFTBufferBase<RankFourTensor> & elastic_buffer);
  /**
   * Helper to populate initial strain buffer
   */
  void populateEpsilonBuffer(FFTBufferBase<RankTwoTensor> & epsilon_buffer);

  /**
   * Advance Fourier epsilon to next iteration
   */
  void advanceReciprocalEpsilon(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                                FFTBufferBase<RankTwoTensor> & stress_buffer,
                                const FFTBufferBase<RankFourTensor> & gamma_hat);
  /**
   * Update sigma in real space for this iteration
   */
  void updateRealSigma(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                       FFTBufferBase<RankTwoTensor> & stress_buffer,
                       FFTBufferBase<RankFourTensor> & elastic_tensor,
                       RankFourTensor & elastic_tensor_homo);

  void filloutElasticTensor(const FFTBufferBase<Real> & ratio_buffer,
                            FFTBufferBase<Real> & index_buffer,
                            FFTBufferBase<RankFourTensor> & elastic_tensor_buffer);
};
