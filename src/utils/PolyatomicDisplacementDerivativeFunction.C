/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "PolyatomicDisplacementDerivativeFunction.h"

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

PolyatomicDisplacementDerivativeFunction::PolyatomicDisplacementDerivativeFunction(
    std::vector<MyTRIM_NS::Element> polyatomic_material,
    nrt_type damage_function_type,
    const PolyatomicDisplacementFunction * net_disp,
    std::vector<std::vector<Real>> Ecap)
  : PolyatomicDisplacementFunctionBase(polyatomic_material, damage_function_type, Ecap),
    _net_displacement_function(net_disp)
{
  if (damage_function_type != NET_DERIVATIVE)
    throw std::exception();

  // set up the gsl ODE machinery
  auto func = &PolyatomicDisplacementDerivativeFunction::odeRHS;
  _sys = {func, NULL, _problem_size, this};
  _ode_driver = gsl_odeiv2_driver_alloc_y_new(&_sys, gsl_odeiv2_step_rk4, 10.0, 1e-2, 1.0e-3);

  Real Edisp_min = std::numeric_limits<Real>::max();
  for (unsigned int j = 0; j < _n_species; ++j)
    if (_material->_element[j]._Edisp < Edisp_min)
      Edisp_min = _material->_element[j]._Edisp;
  _energy_history[0] = Edisp_min;

  // note that initial conditions are 0s for theta_ijl
}

int
PolyatomicDisplacementDerivativeFunction::odeRHS(Real energy,
                                                 const Real disp[],
                                                 Real f[],
                                                 void * params)
{
  (void)disp;
  PolyatomicDisplacementDerivativeFunction * padf =
      (PolyatomicDisplacementDerivativeFunction *)params;
  for (unsigned int i = 0; i < padf->nSpecies(); ++i)
  {
    Real stopping_power = padf->stoppingPower(i, energy);
    for (unsigned int j = 0; j < padf->nSpecies(); ++j)
      for (unsigned int l = 0; l < padf->nSpecies(); ++l)
      {
        // working on the right hand side for theta_ijl
        unsigned int n = padf->mapIndex(i, j, l);
        f[n] = 0;

        for (unsigned int k = 0; k < padf->nSpecies(); ++k)
          f[n] +=
              padf->numberFraction(k) *
              (padf->integralTypeI(energy, i, j, l, k) + padf->integralTypeII(energy, i, j, l, k)) /
              stopping_power;
        f[n] += padf->source(energy, i, j, l) / stopping_power;
      }
  }
  return GSL_SUCCESS;
}

Real
PolyatomicDisplacementDerivativeFunction::integralTypeI(
    Real energy, unsigned int i, unsigned int j, unsigned int l, unsigned int k) const
{
  Real upper_integration_limit = energy * _lambda[i][k];
  Real integral = 0;
  Real Ebind = _material->_element[k]._Elbind;

  // the integration follows the already existing energies
  for (unsigned int n = 1; n < _energy_history.size(); ++n)
  {
    // adjust integration limits to be within E_{l-1} ... min(E_l, upper_integration_limit)
    Real lower = _energy_history[n - 1];
    Real upper = std::min(_energy_history[n], upper_integration_limit);

    // nothing to integrate
    if (lower > upper)
      continue;

    // now integrate from lower to upper
    Real f = 0.5 * (upper - lower);
    for (unsigned int qp = 0; qp < _quad_order; ++qp)
    {
      Real recoil_energy = f * (_quad_points[qp] + 1) + lower;
      integral += f * _quad_weights[qp] * scatteringCrossSection(i, k, energy, recoil_energy) *
                  displacementProbability(k, recoil_energy) *
                  linearInterpolation(k, j, l, recoil_energy - Ebind);
    }
  }
  return integral;
}

Real
PolyatomicDisplacementDerivativeFunction::integralTypeII(
    Real energy, unsigned int i, unsigned int j, unsigned int l, unsigned int k) const
{
  Real upper_integration_limit = energy * _lambda[i][k];
  Real threshold = std::min(_asymptotic_threshold, energy * _lambda[i][k]);

  // store the current displacement function value
  Real current_value = linearInterpolation(i, j, l, energy);

  // estimate the derivative d(nu_i) / dE:
  Real dE = _energy_history.back() - _energy_history[_energy_history.size() - 2];
  Real derivative = (current_value - linearInterpolation(i, j, l, energy - dE)) / dE;

  // integrate up to threshold and multiply by estimate of the derivative
  Real integral = -weightedScatteringIntegral(energy, threshold, i, k) * derivative;

  if (energy * _lambda[i][k] <= _asymptotic_threshold)
    return integral;

  // the integration follows the already existing energies
  for (unsigned int n = 0; n < _energy_history.size(); ++n)
  {
    // adjust integration limits to be within
    Real lower;
    if (n == 0)
      lower = _asymptotic_threshold;
    else
      lower = std::max(_energy_history[n - 1], _asymptotic_threshold);

    Real upper = std::min(_energy_history[n], upper_integration_limit);

    // nothing to integrate
    if (lower > upper)
      continue;

    // now integrate from lower to upper
    Real f = 0.5 * (upper - lower);
    for (unsigned int qp = 0; qp < _quad_order; ++qp)
    {
      Real recoil_energy = f * (_quad_points[qp] + 1) + lower;
      integral += f * _quad_weights[qp] * scatteringCrossSection(i, k, energy, recoil_energy) *
                  (nonCaptureProbability(i, k, energy, recoil_energy) *
                       linearInterpolation(i, j, l, energy - recoil_energy) -
                   current_value);
    }
  }
  return integral;
}

Real
PolyatomicDisplacementDerivativeFunction::source(Real energy,
                                                 unsigned int i,
                                                 unsigned int j,
                                                 unsigned int l)
{
  // compute the derivative of the net displacement function at E
  unsigned int index = _net_displacement_function->energyIndex(energy);
  Real d_net_dE = 0;
  // if energy == Emin => index = 0, the net displacement rate is constant below that value
  if (index != 0)
  {
    Real e1 = _net_displacement_function->energyPoint(index - 1);
    Real e2 = _net_displacement_function->energyPoint(index);
    Real v1 = _net_displacement_function->linearInterpolation(i, j, e1);
    Real v2 = _net_displacement_function->linearInterpolation(i, j, e2);
    d_net_dE = (v2 - v1) / (e2 - e1);
  }
  return _net_displacement_function->integralTypeI(energy, i, j, l) +
         _net_displacement_function->integralTypeII(energy, i, j, l) -
         d_net_dE * stoppingPowerDerivative(i, l, energy);
}

Real
PolyatomicDisplacementDerivativeFunction::linearInterpolation(unsigned int i,
                                                              unsigned int j,
                                                              unsigned int l,
                                                              Real energy) const
{
  unsigned int index = energyIndex(energy);
  if (index == 0)
    return _displacement_function[0][mapIndex(i, j, l)];

  return linearInterpolation(i, j, l, energy, index);
}

Real
PolyatomicDisplacementDerivativeFunction::linearInterpolation(
    unsigned int i, unsigned int j, unsigned int l, Real energy, unsigned int index) const
{
  unsigned int k = mapIndex(i, j, l);

  // linear interpolation
  Real e1 = _energy_history[index - 1];
  Real e2 = _energy_history[index];
  Real v1 = _displacement_function[index - 1][k];
  Real v2 = _displacement_function[index][k];

  return v1 + (energy - e1) / (e2 - e1) * (v2 - v1);
}

#endif
