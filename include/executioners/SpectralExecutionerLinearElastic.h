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

  /// Reference Young's modulus
  const Real _young_modulus;

  /// Poisson ratio
  const Real _poisson_ratio;

  /// Current time
  Real _t_current;

  /// Initial strain tensor
  RankTwoTensor _initial_strain_tensor;

  /// Average factor to obtain homogeneous material
  const Real _average_factor;

  /// User-prescribed error for fixed iteration solver
  const Real _solver_error;

private:
  /**
   * Helper function to get the diffusion equation Green's function corresponding to one time step.
   * @param gamma_hat Linear elastic Green operator in Fourier space
   * @param elasticity_tensor Elasticity tensor used to replace Gamma for some frequencies
   */
  void getGreensFunction(FFTBufferBase<RankFourTensor> & gamma_hat,
                         const RankFourTensor & elasticity_tensor);

  /**
   * Compute initial stress based on homogeneous strain and space-varying elasticity tensor.
   * @param epsilon_buffer Small strain FFT buffer
   * @param elastic_buffer Space-varying elasticity tensor FFT buffer
   */
  FFTBufferBase<RankTwoTensor> & getInitialStress(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                                                  FFTBufferBase<RankFourTensor> & elastic_buffer);
  /**
   * Initialize epsilon buffer with homogeneous/global strain.
   * @param epsilon_buffer Small strain FFT buffer
   */
  void fillOutEpsilonBuffer(FFTBufferBase<RankTwoTensor> & epsilon_buffer);

  /**
   * Initialize epsilon buffer with homogeneous/global strain.
   * @param epsilon_buffer Small strain FFT buffer
   */
  void advanceReciprocalEpsilon(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                                FFTBufferBase<RankTwoTensor> & stress_buffer,
                                const FFTBufferBase<RankFourTensor> & gamma_hat);
  /**
   * Update real stress buffer for current iteration.
   * @param epsilon_buffer Small strain FFT buffer
   * @param stress_buffer Stress FFT buffer
   * @param elastic_tensor Small strain FFT buffer
   */
  void updateRealSigma(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                       FFTBufferBase<RankTwoTensor> & stress_buffer,
                       FFTBufferBase<RankFourTensor> & elastic_tensor);

  /**
   * Initialize space-varying elastic coefficient tensor.
   * @param ratio_buffer Initial condition ratio to distribute material's stiffnesses
   * @param index_buffer Real or integer buffer assigning index to each cell
   * @param elastic_tensor_buffer Space-varying elastic coefficient tensor
   */
  void filloutElasticTensor(const FFTBufferBase<Real> & ratio_buffer,
                            FFTBufferBase<Real> & index_buffer,
                            FFTBufferBase<RankFourTensor> & elastic_tensor_buffer);
  /**
   * Convergence check base on global equilibrium.
   * @param stress Stress tensor for convergence check
   */
  bool hasStressConverged(const FFTBufferBase<RankTwoTensor> & stress);
};
