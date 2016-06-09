/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "DiscretePKAPDF.h"
#include "MooseError.h"
#include "MooseRandom.h"
#include "MagpieUtils.h"

DiscretePKAPDF::DiscretePKAPDF(Real magnitude, std::vector<unsigned int> ZAID, std::vector<Real> energies,
                               unsigned int na, unsigned int np, MultiIndex<Real> probabilities) :
    DiscretePKAPDFBase(magnitude),
    _zaids(ZAID),
    _nZA(_zaids.size()),
    _energies(energies),
    _ng(_energies.size() - 1),
    _na(na),
    _dphi(2.0 * libMesh::pi / _na),
    _np(np),
    _dmu(2.0 / _np),
    _probabilities(probabilities),
    _marginal_cdf_mu(probabilities),
    _marginal_cdf_phi(probabilities),
    _marginal_cdf_energy(probabilities),
    _marginal_cdf_zaid(probabilities)
{
  // check the size of _probabilities
  if (_probabilities.dim() != 4)
    mooseError("probabilities MultiIndex object has wrong dimensions.");

  std::vector<unsigned int> shape = _probabilities.size();
  if (shape[0] != _nZA || shape[1] != _ng || shape[2] != _na || shape[3] != _np)
    mooseError("Size of probabilities is inconsistent with random variable input.");

  precomputeCDF();
}

void
DiscretePKAPDF::precomputeCDF()
{
  /**
   * Compute the weighted pdf: pdf value * width of the bin
   */
  std::vector<unsigned int> shape(4);
  std::vector<unsigned int> index(4);
  for (MultiIndex<Real>::iterator it = _probabilities.begin(); it != _probabilities.end(); ++it)
  {
    index = it.indices().first;
    Real wt = _dmu * _dphi * (_energies[index[1] + 1] - _energies[index[1]]);
    *it *= wt;
  }

  /**
   *  precompute _marginal_cdf_mu
   */
  // Step 1: Marginalize the pdf
  shape.assign(1, _nZA);
  _marginal_cdf_zaid = MultiIndex<Real>(shape);
  for (unsigned int jZA = 0; jZA < _nZA; ++jZA)
  {
    std::vector<unsigned int> local_index(1);
    Real integral_value = 0.0;
    for (unsigned int jE = 0; jE < _ng; ++jE)
      for (unsigned int jP = 0; jP < _na; ++jP)
        for (unsigned int jM = 0; jM < _np; ++jM)
        {
          index[0] = jZA;
          index[1] = jE;
          index[2] = jP;
          index[3] = jM;
          integral_value += _probabilities(index);
        }
    local_index[0] = jZA;
    _marginal_cdf_zaid(local_index) = integral_value;
  }

  // Step 2: Compute cdf
  for (MultiIndex<Real>::iterator it_zaid = _marginal_cdf_zaid.begin(); it_zaid != _marginal_cdf_zaid.end(); ++it_zaid)
  {
    index = it_zaid.indices().first;
    index[0] -= 1;
    if (it_zaid.indices().first[0] > 0)
      *it_zaid += _marginal_cdf_zaid(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (MultiIndex<Real>::iterator it_zaid = _marginal_cdf_zaid.begin(); it_zaid != _marginal_cdf_zaid.end(); ++it_zaid)
  {
    index = it_zaid.indices().first;
    index[0] = _nZA - 1;
    *it_zaid /= _marginal_cdf_zaid(index);
  }

  /**
   *  precompute _marginal_cdf_energy
   */
  // Step 1: Marginalize the pdf
  shape.resize(2);
  index.resize(4);
  shape[0] = _nZA;
  shape[1] = _ng;
  _marginal_cdf_energy = MultiIndex<Real>(shape);
  for (unsigned int jZA = 0; jZA < _nZA; ++jZA)
    for (unsigned int jE = 0; jE < _ng; ++jE)
    {
      std::vector<unsigned int> local_index(2);
      Real integral_value = 0.0;
      for (unsigned int jP = 0; jP < _na; ++jP)
        for (unsigned int jM = 0; jM < _np; ++jM)
        {
          index[0] = jZA;
          index[1] = jE;
          index[2] = jP;
          index[3] = jM;
          integral_value += _probabilities(index);
        }
      local_index[0] = jZA;
      local_index[1] = jE;
      _marginal_cdf_energy(local_index) = integral_value;
    }

  // Step 2: Compute cdf
  for (MultiIndex<Real>::iterator it_energy = _marginal_cdf_energy.begin(); it_energy != _marginal_cdf_energy.end(); ++it_energy)
  {
    index = it_energy.indices().first;
    index[1] -= 1;
    if (it_energy.indices().first[1] > 0)
      *it_energy += _marginal_cdf_energy(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (MultiIndex<Real>::iterator it_energy = _marginal_cdf_energy.begin(); it_energy != _marginal_cdf_energy.end(); ++it_energy)
  {
    index = it_energy.indices().first;
    index[1] = _ng - 1;
    *it_energy /= _marginal_cdf_energy(index);
  }

  /**
   *  precompute _marginal_cdf_phi
   */
  // Step 1: Marginalize the pdf
  shape.resize(3);
  index.resize(4);
  shape[0] = _nZA;
  shape[1] = _ng;
  shape[2] = _na;
  _marginal_cdf_phi = MultiIndex<Real>(shape);
  for (unsigned int jZA = 0; jZA < _nZA; ++jZA)
    for (unsigned int jE = 0; jE < _ng; ++jE)
      for (unsigned int jP = 0; jP < _na; ++jP)
      {
        std::vector<unsigned int> local_index(3);
        Real integral_value = 0.0;
        for (unsigned int jM = 0; jM < _np; ++jM)
        {
          index[0] = jZA;
          index[1] = jE;
          index[2] = jP;
          index[3] = jM;
          integral_value += _probabilities(index);
        }
        local_index[0] = jZA;
        local_index[1] = jE;
        local_index[2] = jP;
        _marginal_cdf_phi(local_index) = integral_value;
      }

  // Step 2: Compute cdf
  for (MultiIndex<Real>::iterator it_phi = _marginal_cdf_phi.begin(); it_phi != _marginal_cdf_phi.end(); ++it_phi)
  {
    index = it_phi.indices().first;
    index[2] -= 1;
    if (it_phi.indices().first[2] > 0)
      *it_phi += _marginal_cdf_phi(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (MultiIndex<Real>::iterator it_phi = _marginal_cdf_phi.begin(); it_phi != _marginal_cdf_phi.end(); ++it_phi)
  {
    index = it_phi.indices().first;
    index[2] = _na - 1;
    *it_phi /= _marginal_cdf_phi(index);
  }

  /**
   *  precompute _marginal_cdf_mu
   */
  // Step 1: No need to marginalize pdf for mu
  _marginal_cdf_mu = _probabilities;

  // Step 2: Compute cdf
  for (MultiIndex<Real>::iterator it_mu = _marginal_cdf_mu.begin(); it_mu != _marginal_cdf_mu.end(); ++it_mu)
  {
    index = it_mu.indices().first;
    index[3] -= 1;
    if (it_mu.indices().first[3] > 0)
      *it_mu += _marginal_cdf_mu(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (MultiIndex<Real>::iterator it_mu = _marginal_cdf_mu.begin(); it_mu != _marginal_cdf_mu.end(); ++it_mu)
  {
    index = it_mu.indices().first;
    index[3] = _np - 1;
    *it_mu /= _marginal_cdf_mu(index);
  }
}

unsigned int
DiscretePKAPDF::sampleHelper(const MultiIndex<Real> & marginal_pdf, const std::vector<unsigned int> indices) const
{
  if (marginal_pdf.dim() - indices.size() != 1)
    mooseError("For sampling the indices vector must reduce the marginal_pdf to a one-dimensional ladder function.");
  std::vector<unsigned int> dimension(indices.size());
  for (unsigned int j = 0; j < indices.size(); ++j)
    dimension[j] = j;
  MultiIndex<Real> new_marginal_pdf = marginal_pdf.slice(dimension, indices);
  return sampleHelper(new_marginal_pdf);
}

unsigned int
DiscretePKAPDF::sampleHelper(const MultiIndex<Real> & marginal_pdf) const
{
  Real r = MooseRandom::rand();
  std::vector<unsigned int> index(1);
  unsigned int j = 0;
  for (; j < marginal_pdf.size()[0]; ++j)
  {
    index[0] = j;
    if (r <= marginal_pdf(index))
      break;
  }
  return j;
}

void
DiscretePKAPDF::drawSample(initialPKAState & initial_state)
{
  /**
   * Same the discrete pdfs by first sampling from the "most marginal"
   * _marginal_cdf_zaid. Then get the conditional _marginal_cdf_energy(zaid)
   * and sample the energy bin from it, then continue on sampling phi...
   */
  std::vector<unsigned int> sampled_indices;
  sampled_indices.push_back(sampleHelper(_marginal_cdf_zaid));
  sampled_indices.push_back(sampleHelper(_marginal_cdf_energy, sampled_indices));
  sampled_indices.push_back(sampleHelper(_marginal_cdf_phi, sampled_indices));
  sampled_indices.push_back(sampleHelper(_marginal_cdf_mu, sampled_indices));

  // first we need to compute Z and m from ZAID
  initial_state._Z = MagpieUtils::getZFromZAID(_zaids[sampled_indices[0]]);
  initial_state._mass = MagpieUtils::getAFromZAID(_zaids[sampled_indices[0]]);

  // the real random variables also need to be resampled uniformly
  // within bin index[j]
  initial_state._energy = (_energies[sampled_indices[1] + 1] - _energies[sampled_indices[1]]) * MooseRandom::rand() + _energies[sampled_indices[1]];

  Real sampled_phi = _dphi * MooseRandom::rand() + _dphi * sampled_indices[2];
  Real sampled_mu = _dmu * MooseRandom::rand() + _dmu * sampled_indices[3] - 1.0;
  // NOTE: mu is measured w.r.t. the x axis because of the convention in Rattlesnake
  // comparing to D&H Eq. 2-41: x -> y, y -> z, z -> x
  initial_state._direction(0) = sampled_mu;
  initial_state._direction(1) = std::sqrt(1.0 - sampled_mu * sampled_mu) * std::cos(sampled_phi);
  initial_state._direction(2) = std::sqrt(1.0 - sampled_mu * sampled_mu) * std::sin(sampled_phi);
}
