/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "SpectralExecutionerLinearElastic.h"
#include "InputParameters.h"
#include "RankFourTensor.h"

registerMooseObject("MagpieApp", SpectralExecutionerLinearElastic);

Real
deltaij(const int i, const int j)
{
  if (i == j)
    return 1.0;
  else
    return 0.0;
}

InputParameters
SpectralExecutionerLinearElastic::validParams()
{
  InputParameters params = SpectralExecutionerBase::validParams();
  params.addClassDescription("Executioner for spectral solve of diffusion equation");
  params.addParam<Real>("time_step", 1.0, "Time step for ODE integration");
  params.addParam<unsigned int>("number_steps", 0.0, "Time step for ODE integration");
  params.addParam<Real>("young_modulus", 1.0, "First parameter for isotropic materials");
  params.addParam<Real>("poisson_ratio", 1.0, "Second parameter for isotropic materials");
  params.addParam<Real>("average_material_factor", 1.0, "Homogeneized factor for multiphase");
  params.addParam<Real>(
      "initial_shear_strain", 0.001, "Homogeneous two-dimensional shear deformation");

  return params;
}

SpectralExecutionerLinearElastic::SpectralExecutionerLinearElastic(
    const InputParameters & parameters)
  : SpectralExecutionerBase(parameters),
    _dt(getParam<Real>("time_step")),
    _nsteps(getParam<unsigned int>("number_steps")),
    _young_modulus(getParam<Real>("young_modulus")),
    _poisson_ratio(getParam<Real>("poisson_ratio")),
    _initial_shear_strain(getParam<Real>("initial_shear_strain")),
    _average_factor(getParam<Real>("average_material_factor"))
{
  // Add check that's a 2D problem with LIBMESH_DIM == 2
  _initial_strain_tensor(0, 0) = 0.0;
  _initial_strain_tensor(0, 1) = _initial_strain_tensor(1, 0) = _initial_shear_strain;
  _initial_strain_tensor(2, 1) = _initial_strain_tensor(1, 2) = 0.0;
  _initial_strain_tensor(0, 2) = _initial_strain_tensor(2, 0) = 0.0;
  _initial_strain_tensor(1, 1) = _initial_strain_tensor(2, 2) = 0.0;

  _t_current = 0.0;
}

void
SpectralExecutionerLinearElastic::populateEpsilonBuffer(
    FFTBufferBase<RankTwoTensor> & epsilon_buffer)
{
  const auto & grid = epsilon_buffer.grid();
  const int ni = grid[0];
  const int nj = grid[1]; // Real space.
  const int nk = grid[2]; // Real space.

  std::size_t index = 0;

  for (int freq_x = 0; freq_x < ni; ++freq_x)
    for (int freq_y = 0; freq_y < nj; ++freq_y)
      for (int freq_z = 0; freq_z < nk; ++freq_z)
      {
        epsilon_buffer.realSpace()[index] = _initial_strain_tensor;
        index++;
      }
}

void
SpectralExecutionerLinearElastic::getGreensFunction(FFTBufferBase<RankFourTensor> & gamma_hat,
                                                    FFTBufferBase<Real> & ratio_buffer,
                                                    const RankFourTensor & elasticity_tensor)
{
  const auto & grid = gamma_hat.grid();
  const int ndim = 3;

  const int ni = grid[0];
  const int nj = grid[1];            // Real space.
  const int nk = (grid[2] >> 1) + 1; // Real space.

  std::size_t index = 0;

  // To plug the right frequencies, they are not right
  const auto & ivec = gamma_hat.kTable(0);
  const auto & jvec = gamma_hat.kTable(1);
  const auto & kvec = gamma_hat.kTable(2);
  /// return the size of the box
  const Point box_size = gamma_hat.getBoxSize();

  //    std::vector<Real> freqs;
  //    for (Real value = -grid[0]/2.0; value < grid[0]/2.0; value+=grid[0]/box_size(0))
  //      freqs.push_back(value);

  /// return the number of grid cells along each dimension without padding
  const std::vector<int> grid_vector = gamma_hat.grid();

  const Complex I(0.0, 1.0);

  for (int freq_x = 0; freq_x < ni; ++freq_x)
    for (int freq_y = 0; freq_y < nj; ++freq_y)
      for (int freq_z = 0; freq_z < nk; ++freq_z)
      {
        const std::array<Complex, 3> freq{ivec[freq_x] * I, jvec[freq_y] * I, kvec[freq_z] * I};

        Real lambda0 = _young_modulus * _average_factor * _poisson_ratio /
                       ((1 + _poisson_ratio) * (1 - 2 * _poisson_ratio));
        Real nu0 = _young_modulus * _average_factor / (2 * (1 + _poisson_ratio));
        Real constant = (lambda0 + nu0) / (nu0 * (lambda0 + 2.0 * nu0));

        for (int i = 0; i < ndim; i++)
          for (int j = 0; j < ndim; j++)
            for (int k = 0; k < ndim; k++)
              for (int l = 0; l < ndim; l++)
              {
                Complex q_square = freq[0] * freq[0] + freq[1] * freq[1] + freq[2] * freq[2];
                if (std::abs(q_square) > 1.0e-12)
                {
                  gamma_hat.reciprocalSpace()[index](i, j, k, l) =
                      -1.0 * constant * (freq[i] * freq[j] * freq[k] * freq[l]) /
                          (q_square * q_square) +
                      (deltaij(i, k) * freq[j] * freq[l] + deltaij(j, k) * freq[i] * freq[l] +
                       deltaij(i, l) * freq[j] * freq[k] + deltaij(j, l) * freq[i] * freq[k]) /
                          (4.0 * nu0 * q_square);
                }
                else
                {
                  gamma_hat.reciprocalSpace()[index] = elasticity_tensor.invSymm();
                }
              }

        index++;
      }
  // exit(1);
}

FFTBufferBase<RankTwoTensor> &
SpectralExecutionerLinearElastic::getInitialStress(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                                                   FFTBufferBase<RankFourTensor> & elastic_tensor)
{
  auto & stress_buffer = getFFTBuffer<RankTwoTensor>("stress");
  const auto & grid = epsilon_buffer.grid();
  int ni = grid[0];
  int nj = grid[1];
  int nk = grid[2];

  // Homogeneous initial state of strain
  populateEpsilonBuffer(epsilon_buffer);

  size_t index = 0;
  for (int freq_x = 0; freq_x < ni; ++freq_x)
    for (int freq_y = 0; freq_y < nj; ++freq_y)
      for (int freq_z = 0; freq_z < nk; ++freq_z)
      {
        stress_buffer.realSpace()[index] =
            elastic_tensor.realSpace()[index] * epsilon_buffer.realSpace()[index];
        index++;
      }

  return stress_buffer;
}

void
SpectralExecutionerLinearElastic::advanceReciprocalEpsilon(
    FFTBufferBase<RankTwoTensor> & epsilon_buffer,
    FFTBufferBase<RankTwoTensor> & stress_buffer,
    const FFTBufferBase<RankFourTensor> & gamma_hat)
{
  const auto & grid = epsilon_buffer.grid();
  int ni = grid[0];
  int nj = grid[1];
  int nk = (grid[2] >> 1) + 1;

  size_t index = 0;
  // To plug the right frequencies, they are not right
  const auto & ivec = gamma_hat.kTable(0);
  const auto & jvec = gamma_hat.kTable(1);
  const auto & kvec = gamma_hat.kTable(2);

  for (int freq_x = 0; freq_x < ni; ++freq_x)
    for (int freq_y = 0; freq_y < nj; ++freq_y)
      for (int freq_z = 0; freq_z < nk; ++freq_z)
      {
        if (ivec[freq_x] == 0 && jvec[freq_y] == 0 && kvec[freq_z] == 0)
        {
          epsilon_buffer.reciprocalSpace()[index] = _initial_strain_tensor;
        }
        else
        {
          epsilon_buffer.reciprocalSpace()[index] =
              gamma_hat.reciprocalSpace()[index] * stress_buffer.reciprocalSpace()[index];
        }
        index++;
      }
}

void
SpectralExecutionerLinearElastic::updateRealSigma(FFTBufferBase<RankTwoTensor> & epsilon_buffer,
                                                  FFTBufferBase<RankTwoTensor> & stress_buffer,
                                                  FFTBufferBase<RankFourTensor> & elastic_tensor,
                                                  RankFourTensor & elastic_tensor_homo)
{
  const auto & grid = epsilon_buffer.grid();
  int ni = grid[0];
  int nj = grid[1];
  int nk = grid[2];

  size_t index = 0;

  for (int freq_x = 0; freq_x < ni; ++freq_x)
    for (int freq_y = 0; freq_y < nj; ++freq_y)
      for (int freq_z = 0; freq_z < nk; ++freq_z)
      {
        stress_buffer.realSpace()[index] =
            (elastic_tensor.realSpace()[index] - elastic_tensor_homo) * epsilon_buffer.realSpace()[index];
        index++;
      }
}

void
SpectralExecutionerLinearElastic::filloutElasticTensor(
    const FFTBufferBase<Real> & ratio_buffer,
    FFTBufferBase<Real> & index_buffer,
    FFTBufferBase<RankFourTensor> & elastic_tensor_buffer)
{
  const auto & grid = elastic_tensor_buffer.grid();

  int ni = grid[0];
  int nj = grid[1];
  int nk = grid[2];

  size_t index = 0;

  for (int freq_x = 0; freq_x < ni; ++freq_x)
    for (int freq_y = 0; freq_y < nj; ++freq_y)
      for (int freq_z = 0; freq_z < nk; ++freq_z)
      {
        // Define elastic tensor
        std::vector<Real> iso_const(2);
        iso_const[0] = _young_modulus * ratio_buffer.realSpace()[index];
        iso_const[1] = _poisson_ratio;

        elastic_tensor_buffer.realSpace()[index].fillFromInputVector(
            iso_const, RankFourTensor::symmetric_isotropic_E_nu);
        index_buffer.realSpace()[index] = index;
        index++;
      }
}

void
SpectralExecutionerLinearElastic::execute()
{
  unsigned int thisStep = 0;

  _time_step = thisStep;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _fe_problem.advanceState();

  auto & epsilon_buffer = getFFTBuffer<RankTwoTensor>("epsilon");

  auto & ratio_buffer = getFFTBuffer<Real>("stiffness_ratio");
  auto & elastic_tensor_buffer = getFFTBuffer<RankFourTensor>("elastic");
  auto & index_buffer = getFFTBuffer<Real>("index_buffer");

  filloutElasticTensor(ratio_buffer, index_buffer, elastic_tensor_buffer);

  // Get corresponding initial stress (also fill out epsilon_buffer with initial strain)
  auto & stress_buffer = getInitialStress(epsilon_buffer, elastic_tensor_buffer);

  // Get specific Green's function
  auto & gamma_hat = getFFTBuffer<RankFourTensor>("gamma");

  if (gamma_hat.dim() != 3)
    mooseError("Error: Dimension not implemented in SpectralExecutionerLinearElastic.");

  thisStep++;
  _t_current += _dt;
  _time_step = thisStep;

  _fe_problem.execute(EXEC_FINAL);
  _time = _t_current;

  Moose::out << "_t_current: " << _t_current << ". \n";

  _fe_problem.outputStep(EXEC_FINAL);

  ///
  std::vector<Real> iso_const(2);
  iso_const[0] = _young_modulus * _average_factor;
  iso_const[1] = _poisson_ratio;
  RankFourTensor elasticity_homo;
  elasticity_homo.fillFromInputVector(iso_const, RankFourTensor::symmetric_isotropic_E_nu);
  getGreensFunction(gamma_hat, ratio_buffer, elasticity_homo);
  ///

  // Our plans do not preserve the inputs (unfortunately)
  FFTData<RankTwoTensor> stress_buffer_backup_real = stress_buffer.realSpace();
  epsilon_buffer.forward();

  FFTData<ComplexRankTwoTensor> epsilon_buffer_backup_reciprocal =
      epsilon_buffer.reciprocalSpace();

  for (unsigned int step_no = 0; step_no < _nsteps; step_no++)
  {
    // Preserve data
    updateRealSigma(epsilon_buffer, stress_buffer, elastic_tensor_buffer, elasticity_homo);
    // We would need here the stress for convergence check
    // (a) plus bookkeeping
    stress_buffer_backup_real = stress_buffer.realSpace();
    stress_buffer.forwardRaw();
    stress_buffer.reciprocalSpace() *= stress_buffer.backwardScale();
    stress_buffer.realSpace() = stress_buffer_backup_real;
    // Compute new strain tensor in Fourier space
    // (c)
    // Preserve data

    epsilon_buffer.reciprocalSpace() = epsilon_buffer_backup_reciprocal;
    advanceReciprocalEpsilon(epsilon_buffer, stress_buffer, gamma_hat);

    // (d)
    epsilon_buffer_backup_reciprocal = epsilon_buffer.reciprocalSpace();
    // For output purposes
    epsilon_buffer.backwardRaw();

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

