/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "PolyatomicDisplacementFunctionBase.h"

// mytrim includes
#include <mytrim/simconf.h>
#include <mytrim/ion.h>
#include <mytrim/element.h>

// general includes
#include <assert.h>
#include <limits>

// gsl includes
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

PolyatomicDisplacementFunctionBase::PolyatomicDisplacementFunctionBase(
    std::vector<MyTRIM_NS::Element> polyatomic_material,
    nrt_type damage_function_type,
    std::vector<std::vector<Real>> Ecap)
  : _damage_function_type(damage_function_type),
    _quad_order(4),
    _n_species(polyatomic_material.size()),
    _problem_size(_damage_function_type == ENERGY ? _n_species
                                                  : _damage_function_type == NET_DERIVATIVE
                                                        ? _n_species * _n_species * _n_species
                                                        : _n_species * _n_species),
    _energy_history({0.0}),
    _displacement_function({std::vector<Real>(_problem_size)})
{
  _simconf = libmesh_make_unique<MyTRIM_NS::SimconfType>(1);
  _material = libmesh_make_unique<MyTRIM_NS::MaterialBase>(_simconf.get(), 1);
  for (auto & elem : polyatomic_material)
  {
    _material->_element.push_back(elem);
    _ions.push_back(libmesh_make_unique<MyTRIM_NS::IonBase>(elem._Z, elem._m, 0));
  }
  _material->prepare();

  // compute _lambda
  _lambda.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    _lambda[i].resize(_n_species);
    for (unsigned int j = 0; j < _n_species; ++j)
    {
      Real Ai = _material->_element[i]._m;
      Real Aj = _material->_element[j]._m;
      _lambda[i][j] = 4.0 * Ai * Aj / (Ai + Aj) / (Ai + Aj);
    }
  }

  if (Ecap.size() == 1 && Ecap[0].size() == 0)
  {
    // default is Ecap_ij = Ed_j
    _Ecap.resize(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      _Ecap[i].resize(_n_species);
      for (unsigned int j = 0; j < _n_species; ++j)
        _Ecap[i][j] = _material->_element[j]._Edisp;
    }
  }
  else
    _Ecap = Ecap;

  // set up integration rule
  auto * qp_table = gsl_integration_glfixed_table_alloc(_quad_order);
  _quad_points.resize(_quad_order);
  _quad_weights.resize(_quad_order);
  for (std::size_t j = 0; j < _quad_order; ++j)
  {
    double point, weight;
    gsl_integration_glfixed_point(-1.0, 1.0, j, &point, &weight, qp_table);
    _quad_points[j] = point;
    _quad_weights[j] = weight;
  }
  gsl_integration_glfixed_table_free(qp_table);
}

PolyatomicDisplacementFunctionBase::~PolyatomicDisplacementFunctionBase()
{
  gsl_odeiv2_driver_free(_ode_driver);
}

void
PolyatomicDisplacementFunctionBase::advanceDisplacements(Real new_energy)
{
  if (_energy_history.back() > new_energy)
    return;

  // add the new energy point to energy history and store old_energy
  Real old_energy = _energy_history.back();
  _energy_history.push_back(new_energy);

  // add new energy step to _displacement_function; initial guess is the old energy
  _displacement_function.push_back(_displacement_function.back());

  int status = gsl_odeiv2_driver_apply(
      _ode_driver, &old_energy, new_energy, _displacement_function.back().data());

  if (status != GSL_SUCCESS)
    std ::cout << "gsl_odeiv2_driver_apply returned error code  " << status << std::endl;
}

void
PolyatomicDisplacementFunctionBase::computeDisplacementFunctionIntegral()
{
  // clear the _displacement_function_integral array because this might
  // have been called before; not the most efficient if the user keeps on calling
  // this function but they are at their own peril and should know to call it
  // after finishing the computation of disp function
  _displacement_function_integral.clear();

  // resize the array; note that initial value is zero (entry for
  // _displacement_function_integral[0])
  _displacement_function_integral.resize(nEnergySteps(), std::vector<Real>(_problem_size));

  for (unsigned int e = 1; e < nEnergySteps(); ++e)
  {
    _displacement_function_integral.push_back(std::vector<Real>(_problem_size));

    // set energy points for integration
    Real lower = energyPoint(e - 1);
    Real upper = energyPoint(e);
    Real f = 0.5 * (upper - lower);
    for (unsigned int qp = 0; qp < _quad_order; ++qp)
    {
      Real energy = f * (_quad_points[qp] + 1) + lower;
      Real w = f * _quad_weights[qp];

      for (unsigned int n = 0; n < _problem_size; ++n)
      {
        unsigned int i, j, l;
        inverseMapIndex(n, i, j, l);
        _displacement_function_integral[e][n] += w * linearInterpolation(energy, i, j, l);
      }
    }
  }
}

Real
PolyatomicDisplacementFunctionBase::stoppingPower(unsigned int species, Real energy)
{
  assert(species < _n_species);
  _ions[species]->_E = energy;
  // need to divide by atomic density because we need Se for PC equations
  // units are eV * Angstrom**2
  return _material->getrstop(_ions[species].get()) / _material->_arho;
}

Real
PolyatomicDisplacementFunctionBase::stoppingPowerDerivative(unsigned int species,
                                                            unsigned int changing_species,
                                                            Real energy)
{
  _ions[species]->_E = energy;
  return _material->getDrstopDcomp(_ions[species].get(), _material->_element[changing_species]);
}

/**
 * The scattering cross section is computed using Lindhard's universal formula
 * with f(xi) evaluated using the Thomas Fermi potential. Units are Angstrom**2 / eV
 */
Real
PolyatomicDisplacementFunctionBase::scatteringCrossSection(unsigned int i,
                                                           unsigned int j,
                                                           Real energy,
                                                           Real recoil_energy) const
{
  assert(i < _n_species);
  assert(j < _n_species);

  // get the current A & Z for projectile: i and target: j
  Real Ai = _material->_element[i]._m;
  Real Aj = _material->_element[j]._m;
  Real Zi = _material->_element[i]._Z;
  Real Zj = _material->_element[j]._Z;

  // Z only appears as 1/3 power, so it's better to apply the 1/3 power here to get a sqrt
  Real Z = std::sqrt(std::pow(Zi, 2.0 / 3.0) + std::pow(Zj, 2.0 / 3.0));
  Real a = 0.8853 * _abohr / Z;

  // compute El and Tm
  Real Tm = _lambda[i][j] * energy;
  Real El = 30.7514664 * Zi * Zj * Z * (Ai + Aj) / Aj;

  // compute sqrt(t) and then the cross section using Lindhard's formula
  Real sqrt_t = std::sqrt((energy * energy / El / El) * recoil_energy / Tm);
  return 0.5 * universalF(sqrt_t) / recoil_energy / sqrt_t * M_PI * a * a;
}

Real
PolyatomicDisplacementFunctionBase::universalF(Real xi) const
{
  Real lp = 1.309;
  return lp * std::pow(xi, 1.0 / 3.0) *
         std::pow(1.0 + std::pow(2.0 * lp * std::pow(xi, 4.0 / 3.0), 2.0 / 3.0), -1.5);
}

unsigned int
PolyatomicDisplacementFunctionBase::findSpeciesIndex(unsigned int atomic_number,
                                                     Real mass_number) const
{
  for (unsigned int j = 0; j < _n_species; ++j)
    if (atomic_number == _material->_element[j]._Z &&
        std::abs(mass_number - _material->_element[j]._m) < 1.0e-10)
      return j;
  return libMesh::invalid_uint;
}

Real
PolyatomicDisplacementFunctionBase::displacementProbability(unsigned int k, Real energy) const
{
  return energy < _material->_element[k]._Edisp ? 0 : 1;
}

Real
PolyatomicDisplacementFunctionBase::nonCaptureProbability(unsigned int i,
                                                          unsigned int k,
                                                          Real energy,
                                                          Real recoil_energy) const
{
  // Parkin & Coulter, JNM 101, (1981)
  return 1 -
         displacementProbability(k, recoil_energy) * (energy - recoil_energy < _Ecap[i][k] ? 1 : 0);
}

unsigned int
PolyatomicDisplacementFunctionBase::energyIndex(Real energy) const
{
  return energy >= _energy_history.back()
             ? nEnergySteps() - 1
             : std::distance(
                   _energy_history.begin(),
                   std::upper_bound(_energy_history.begin(), _energy_history.end(), energy));
}

/**
 * weightedScatteringIntegral computes the integral:
 * int_0^{energy_limit} d(sigma_ij) / dT T dT
 */
Real
PolyatomicDisplacementFunctionBase::weightedScatteringIntegral(Real energy,
                                                               Real energy_limit,
                                                               unsigned int i,
                                                               unsigned int j) const
{
  // uses the t -> 0 limit to WSS form of Lindhard's universal cross section with Thomas-Fermi
  // potential for recoil_energy < _asymptotic_threshold
  Real limit = std::min(_asymptotic_threshold, energy_limit);
  Real Mi = _material->_element[i]._m;
  Real Mj = _material->_element[j]._m;
  Real Zi = _material->_element[i]._Z;
  Real Zj = _material->_element[j]._Z;
  Real Z = std::pow(std::pow(Zi, 2.0 / 3.0) + std::pow(Zj, 2.0 / 3.0), 1.5);
  Real a = _abohr * 0.8853 * std::pow(Z, -1.0 / 3.0);
  Real alpha_ij = 0.5 * 1.309 * M_PI * a * a * std::pow(Mi / Mj, 1.0 / 3.0) *
                  std::pow(Zi * Zj * 61.5038 * std::pow(Z, 1.0 / 3.0), 2.0 / 3.0);
  Real integral = 1.5 * std::pow(limit, 2.0 / 3.0) * alpha_ij * std::pow(energy, -1.0 / 3.0);

  if (energy_limit <= _asymptotic_threshold)
    return integral;

  // numerical integration from _asymptotic_threshold to energy_limit
  Real spacing = std::max(energy_limit / _asymptotic_threshold /
                              std::ceil(energy_limit / _asymptotic_threshold / 1.05),
                          1.02);

  std::vector<Real> energies = {_asymptotic_threshold};
  Real current_energy = _asymptotic_threshold;
  for (;;)
  {
    current_energy *= spacing;
    if (current_energy >= energy_limit)
    {
      energies.push_back(energy_limit);
      break;
    }
    energies.push_back(current_energy);
  }

  for (unsigned int l = 1; l < energies.size(); ++l)
  {
    Real lower = energies[l - 1];
    Real upper = energies[l];
    Real scale = 0.5 * (upper - lower);
    for (unsigned int qp = 0; qp < _quad_order; ++qp)
    {
      Real recoil_energy = scale * (_quad_points[qp] + 1) + lower;
      integral += scale * _quad_weights[qp] * scatteringCrossSection(i, j, energy, recoil_energy) *
                  recoil_energy;
    }
  }
  return integral;
}

Real
PolyatomicDisplacementFunctionBase::linearInterpolation(Real energy,
                                                        unsigned int i,
                                                        unsigned int j,
                                                        unsigned int l) const
{
  unsigned int index = energyIndex(energy);
  if (index == 0)
    return _displacement_function[0][mapIndex(i, j, l)];

  return linearInterpolationHelper(energy, index, i, j, l);
}

Real
PolyatomicDisplacementFunctionBase::linearInterpolationHelper(
    Real energy, unsigned int index, unsigned int i, unsigned int j, unsigned int l) const
{
  unsigned int k = mapIndex(i, j, l);

  // linear interpolation
  Real e1 = _energy_history[index - 1];
  Real e2 = _energy_history[index];
  Real v1 = _displacement_function[index - 1][k];
  Real v2 = _displacement_function[index][k];

  return v1 + (energy - e1) / (e2 - e1) * (v2 - v1);
}

Real
PolyatomicDisplacementFunctionBase::linearInterpolationIntegralDamageFunction(Real energy,
                                                                              unsigned int i,
                                                                              unsigned int j,
                                                                              unsigned int l) const
{
  unsigned int index = energyIndex(energy);
  unsigned int k = mapIndex(i, j, l);

  if (index == 0)
    return _displacement_function_integral[0][mapIndex(i, j, l)];

  // linear interpolation
  Real e1 = _energy_history[index - 1];
  Real e2 = _energy_history[index];
  Real v1 = _displacement_function_integral[index - 1][k];
  Real v2 = _displacement_function_integral[index][k];

  return v1 + (energy - e1) / (e2 - e1) * (v2 - v1);
}

#endif
