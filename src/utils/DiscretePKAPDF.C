#include "DiscretePKAPDF.h"
#include "MooseError.h"
#include "MooseRandom.h"
#include "MagpieUtils.h"

DiscretePKAPDF::DiscretePKAPDF(const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities) :
    DiscretePKAPDFBase(ZAID, energies),
    _na(probabilities.size()[2]),
    _dphi(2.0 * libMesh::pi / _na),
    _np(probabilities.size()[3]),
    _dmu(2.0 / _np),
    _marginal_cdf_mu(probabilities),
    _marginal_cdf_phi(probabilities),
    _marginal_cdf_energy(probabilities),
    _marginal_cdf_zaid(probabilities)
{
  // check the size of _probabilities
  if (probabilities.dim() != 4)
    mooseError("probabilities MultiIndex object has wrong dimensions.");

  MultiIndex<Real>::size_type shape = probabilities.size();
  if (shape[0] != _nZA || shape[1] != _ng || shape[2] != _na || shape[3] != _np)
    mooseError("Size of probabilities is inconsistent with random variable input.");

  precomputeCDF(probabilities);
  computeMagnitude(probabilities);
}

void
DiscretePKAPDF::precomputeCDF(MultiIndex<Real> probabilities)
{
  /**
   * Compute the weighted pdf: pdf value * width of the bin
   */
  MultiIndex<Real>::size_type shape(4);
  MultiIndex<Real>::size_type index(4);
  for (auto it: probabilities)
  {
    index = it.first;
    Real wt = _dmu * _dphi * (_energies[index[1] + 1] - _energies[index[1]]);
    it.second *= wt;
  }

  /**
   *  precompute _marginal_cdf_mu
   */
  // Step 1: Marginalize the pdf
  shape.assign(1, _nZA);
  _marginal_cdf_zaid = MultiIndex<Real>(shape);
  for (unsigned int jZA = 0; jZA < _nZA; ++jZA)
  {
    Real integral_value = 0.0;
    for (unsigned int jE = 0; jE < _ng; ++jE)
      for (unsigned int jP = 0; jP < _na; ++jP)
        for (unsigned int jM = 0; jM < _np; ++jM)
          integral_value += probabilities({jZA, jE, jP, jM});

    _marginal_cdf_zaid({jZA}) = integral_value;
  }

  // Step 2: Compute cdf
  for (auto zaid: _marginal_cdf_zaid)
  {
    index = zaid.first;
    index[0] -= 1;
    if (zaid.first[0] > 0)
      zaid.second += _marginal_cdf_zaid(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (auto zaid: _marginal_cdf_zaid)
  {
    index = zaid.first;
    index[0] = _nZA - 1;
    zaid.second /= _marginal_cdf_zaid(index);
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
      Real integral_value = 0.0;
      for (unsigned int jP = 0; jP < _na; ++jP)
        for (unsigned int jM = 0; jM < _np; ++jM)
          integral_value += probabilities({jZA, jE, jP, jM});
      _marginal_cdf_energy({jZA, jE}) = integral_value;
    }

  // Step 2: Compute cdf
  for (auto energy: _marginal_cdf_energy)
  {
    index = energy.first;
    index[1] -= 1;
    if (energy.first[1] > 0)
      energy.second += _marginal_cdf_energy(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (auto energy: _marginal_cdf_energy)
  {
    index = energy.first;
    index[1] = _ng - 1;
    energy.second /= _marginal_cdf_energy(index);
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
        Real integral_value = 0.0;
        for (unsigned int jM = 0; jM < _np; ++jM)
          integral_value += probabilities({jZA, jE, jP, jM});
        _marginal_cdf_phi({jZA, jE, jP}) = integral_value;
      }

  // Step 2: Compute cdf
  for (auto phi: _marginal_cdf_phi)
  {
    index = phi.first;
    index[2] -= 1;
    if (phi.first[2] > 0)
      phi.second += _marginal_cdf_phi(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (auto phi: _marginal_cdf_phi)
  {
    index = phi.first;
    index[2] = _na - 1;
    phi.second /= _marginal_cdf_phi(index);
  }

  /**
   *  precompute _marginal_cdf_mu
   */
  // Step 1: No need to marginalize pdf for mu
  _marginal_cdf_mu = probabilities;

  // Step 2: Compute cdf
  for (auto mu: _marginal_cdf_mu)
  {
    index = mu.first;
    index[3] -= 1;
    if (mu.first[3] > 0)
      mu.second += _marginal_cdf_mu(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (auto mu: _marginal_cdf_mu)
  {
    index = mu.first;
    index[3] = _np - 1;
    mu.second /= _marginal_cdf_mu(index);
  }
}

void
DiscretePKAPDF::drawSample(std::vector<MyTRIM_NS::IonBase> & initial_state) const
{
  //resize initial_state
  initial_state.resize(1);

  /**
   * Same the discrete pdfs by first sampling from the "most marginal"
   * _marginal_cdf_zaid. Then get the conditional _marginal_cdf_energy(zaid)
   * and sample the energy bin from it, then continue on sampling phi...
   */
  MultiIndex<Real>::size_type sampled_indices;
  sampled_indices.push_back(sampleHelper(_marginal_cdf_zaid));
  sampled_indices.push_back(sampleHelper(_marginal_cdf_energy, sampled_indices));
  sampled_indices.push_back(sampleHelper(_marginal_cdf_phi, sampled_indices));
  sampled_indices.push_back(sampleHelper(_marginal_cdf_mu, sampled_indices));

  // first we need to compute Z and m from ZAID
  initial_state[0]._Z = MagpieUtils::getZFromZAID(_zaids[sampled_indices[0]]);
  initial_state[0]._m = MagpieUtils::getAFromZAID(_zaids[sampled_indices[0]]);

  // the real random variables also need to be resampled uniformly
  // within bin index[j]
  initial_state[0]._E = (_energies[sampled_indices[1] + 1] - _energies[sampled_indices[1]]) * MooseRandom::rand() + _energies[sampled_indices[1]];

  // NOTE: we need to sample the direction in this class because the direction is anisotropic
  Real sampled_phi = _dphi * MooseRandom::rand() + _dphi * sampled_indices[2];
  Real sampled_mu = _dmu * MooseRandom::rand() + _dmu * sampled_indices[3] - 1.0;
  initial_state[0]._dir(0) = std::sqrt(1.0 - sampled_mu * sampled_mu) * std::cos(sampled_phi);
  initial_state[0]._dir(1) = std::sqrt(1.0 - sampled_mu * sampled_mu) * std::sin(sampled_phi);
  initial_state[0]._dir(2) = sampled_mu;
}

void
DiscretePKAPDF::computeMagnitude(MultiIndex<Real> probabilities)
{
  _magnitude = 0.0;
  for (auto it : probabilities)
  {
    MultiIndex<Real>::size_type index = it.first;
    Real delE = _energies[index[1] + 1] - _energies[index[1]];
    _magnitude += delE * _dphi * _dmu * it.second;
  }
}
