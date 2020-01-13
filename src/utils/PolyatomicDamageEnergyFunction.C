/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "PolyatomicDamageEnergyFunction.h"

// mytrim includes
#include <mytrim/simconf.h>
#include <mytrim/ion.h>
#include <mytrim/element.h>

// general includes
#include <assert.h>
#include <limits>
#include <exception>

// gsl includes
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

PolyatomicDamageEnergyFunction::PolyatomicDamageEnergyFunction(
    std::vector<MyTRIM_NS::Element> polyatomic_material,
    nrt_type damage_function_type,
    std::vector<std::vector<Real>> Ecap)
  : PolyatomicDisplacementFunctionBase(polyatomic_material, damage_function_type, Ecap)
{
  if (damage_function_type != ENERGY)
    throw std::exception();

  // set the number of indices
  _n_indices = 1;

  /*
   * The ENERGY mode has no lower cutoff, because is nu_i(0) = 0. However, this will
   * lead to issues with evaluating the scattering XS. We use Lindhard's initial condition
   * stated as nu(E) / E -> 1 as E -> 0 [Integral Equations Governing Radiation Effects,
   * Mat. Fys . Medd . Dan . Vid. Selsk . 33, no . 10 (1963)].
   * Threshold is set at 0.01 eV
   */
  _energy_history[0] = 0.01;

  for (unsigned int i = 0; i < _n_species; ++i)
    _displacement_function[0][i] = _energy_history[0];
}

std::vector<Real>
PolyatomicDamageEnergyFunction::getRHS(Real energy)
{
  std::vector<Real> f(_problem_size);
  for (unsigned int i = 0; i < nSpecies(); ++i)
  {
    f[i] = 0;
    Real denominator = stoppingPower(i, energy);

    // range of energies from the threshold to E
    for (unsigned int j = 0; j < nSpecies(); ++j)
    {
      f[i] += numberFraction(j) * integralTypeI(energy, i, j);

      if (energy <= taylorSeriesThreshold())
        denominator +=
            numberFraction(j) * weightedScatteringIntegral(energy, energy * lambda(i, j), i, j);
      else
        f[i] += numberFraction(j) * integralTypeII(energy, i, j);
    }
    f[i] /= denominator;
  }
  return f;
}

Real
PolyatomicDamageEnergyFunction::integralTypeI(Real energy, unsigned int i, unsigned int j)
{
  Real upper_integration_limit = energy * _lambda[i][j];
  Real threshold = std::min(_asymptotic_threshold, energy * _lambda[i][j]);

  // integrate up to threshold, 0, ..., threshold
  Real integral = weightedScatteringIntegral(energy, threshold, i, j);

  if (energy * _lambda[i][j] <= _asymptotic_threshold)
    return integral;

  // start at asymptotic threshold and integrate up to energy
  // the integration follows the already existing energies
  for (unsigned int l = 0; l < _energy_history.size(); ++l)
  {
    Real lower;
    if (l == 0)
      lower = _asymptotic_threshold;
    else
      lower = std::max(_energy_history[l - 1], _asymptotic_threshold);

    Real upper = std::min(_energy_history[l], upper_integration_limit);

    // nothing to integrate
    if (lower > upper)
      continue;

    // now integrate from lower to upper
    Real scale = 0.5 * (upper - lower);
    for (unsigned int qp = 0; qp < _quad_order; ++qp)
    {
      Real recoil_energy = scale * (_quad_points[qp] + 1) + lower;
      integral += scale * _quad_weights[qp] * scatteringCrossSection(i, j, energy, recoil_energy) *
                  linearInterpolation(recoil_energy, j);
    }
  }
  return integral;
}

Real
PolyatomicDamageEnergyFunction::integralTypeII(Real energy, unsigned int i, unsigned int j)
{
  Real upper_integration_limit = energy * _lambda[i][j];
  Real threshold = std::min(_asymptotic_threshold, energy * _lambda[i][j]);

  // store the current displacement function value
  Real current_value = linearInterpolation(energy, i);

  // estimate the derivative d(nu_i) / dE:
  Real dE = _energy_history.back() - _energy_history[_energy_history.size() - 2];
  Real derivative = (current_value - linearInterpolation(energy - dE, i)) / dE;

  // integrate up to threshold and multiply by estimate of the derivative
  Real integral = -weightedScatteringIntegral(energy, threshold, i, j) * derivative;

  if (energy * _lambda[i][j] <= _asymptotic_threshold)
    return integral;

  // integrate from threshold up to energy
  for (unsigned int l = 0; l < _energy_history.size(); ++l)
  {
    Real lower;
    if (l == 0)
      lower = _asymptotic_threshold;
    else
      lower = std::max(_energy_history[l - 1], _asymptotic_threshold);

    Real upper = std::min(_energy_history[l], upper_integration_limit);

    // nothing to integrate
    if (lower > upper)
      continue;

    // now integrate from lower to upper
    Real scale = 0.5 * (upper - lower);
    for (unsigned int qp = 0; qp < _quad_order; ++qp)
    {
      Real recoil_energy = scale * (_quad_points[qp] + 1) + lower;
      integral += scale * _quad_weights[qp] * scatteringCrossSection(i, j, energy, recoil_energy) *
                  (linearInterpolation(energy - recoil_energy, i) - current_value);
    }
  }
  return integral;
}

#endif
