#include "DiscreteFissionPKAPDF.h"
#include "MooseError.h"
#include "MooseRandom.h"
#include "MagpieUtils.h"
#include "MooseUtils.h"
#include <fstream>
#include <map>

DiscreteFissionPKAPDF::DiscreteFissionPKAPDF(Real magnitude, const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities) :
    DiscretePKAPDFBase(magnitude, ZAID, energies),
    _marginal_cdf_target(probabilities),
    _conditional_cdf_energy(probabilities)
{
  // check the size of _probabilities
  if (probabilities.dim() != 2)
    mooseError("probabilities MultiIndex object has wrong dimensions.");

  MultiIndex<Real>::size_type shape = probabilities.size();
  if (shape[0] != _nZA || shape[1] != _ng)
    mooseError("Size of probabilities is inconsistent with random variable input.");

  precomputeCDF(probabilities);

  readFissionData(ZAID);
}

void
DiscreteFissionPKAPDF::precomputeCDF(MultiIndex<Real> probabilities)
{
  // Step 1: Compute the marginal distribution w.r.t. incident energy of the neutron
  MultiIndex<Real>::size_type shape(1);
  MultiIndex<Real>::size_type index(2);

  shape[0] = _nZA;
  _marginal_cdf_target = MultiIndex<Real>(shape);
  for (unsigned int jZA = 0; jZA < _nZA; ++jZA)
  {
    Real integral_value = 0.0;
    for (unsigned int jE = 0; jE < _ng; ++jE)
      integral_value += probabilities({jZA, jE});
    _marginal_cdf_target({jZA}) = integral_value;
  }

  // Step 2: Compute cdf
  for (auto target: _marginal_cdf_target)
  {
    index = target.first;
    index[0] -= 1;
    if (target.first[0] > 0)
      target.second += _marginal_cdf_target(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  index[0] = _nZA - 1;
  Real last_value = _marginal_cdf_target(index);
  for (auto target: _marginal_cdf_target)
    target.second /= last_value;

  // Step 1: Compute the conditional distribution of the target w.r.t the sampled energy
  _conditional_cdf_energy = probabilities;

  // Step 2: Compute cdf
  for (auto energy: _conditional_cdf_energy)
  {
    index = energy.first;
    index[1] -= 1;
    if (energy.first[1] > 0)
      energy.second += _conditional_cdf_energy(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  for (auto energy: _conditional_cdf_energy)
  {
    index = energy.first;
    index[0] = _nZA - 1;
    energy.second /= _conditional_cdf_energy(index);
  }
}

void
DiscreteFissionPKAPDF::drawSample(std::vector<initialPKAState> & initial_state)
{
  //resize initial_state
  initial_state.resize(2);

  // first sample the target ZAID and the energy group of the incident neutron
  std::vector<unsigned int> sampled_indices;
  sampled_indices.push_back(sampleHelper(_marginal_cdf_target));
  sampled_indices.push_back(sampleHelper(_conditional_cdf_energy, sampled_indices));

  // first we need to compute Z and m from ZAID
  auto target_Z = MagpieUtils::getZFromZAID(_zaids[sampled_indices[0]]);
  auto target_A = MagpieUtils::getAFromZAID(_zaids[sampled_indices[0]]);

  // the real random variables also need to be resampled uniformly
  // within bin index[j]
  auto neutron_energy = (_energies[sampled_indices[1] + 1] - _energies[sampled_indices[1]]) * MooseRandom::rand() + _energies[sampled_indices[1]];

  // pull the right fission yield table and sample Z, A, and energy
  auto energy = MagpieUtils::determineNeutronType(neutron_energy);
  auto key = _zaids[sampled_indices[0]];

  // pull the map based on the neutron energy
  auto zaid_map = _fission_zaids[energy];
  auto cdf_map = _fission_cdf[energy];

  /// check if key exists
  if (zaid_map.find(key) == zaid_map.end())
    mooseError("Fission product ZAID's not found for this target isotope, " << key << " at this energy, " << energy);

  if (cdf_map.find(key) == cdf_map.end())
    mooseError("Sum yield not found for this target isotope, " << key << " at this energy, " << energy);

  // now we know the key1/key2 pair exists
  auto zaid_products = zaid_map[key];
  auto cdf_products = cdf_map[key];

  // calculate the total kinetic energy of the fission products
  auto total_ke = determineFragmentsEnergy(target_Z, target_A);

  // sample the fission yield cdf to determine the first fission product
  auto ffindex = sampleHelper(cdf_products);
  auto ffzaid = zaid_products[ffindex];
  initial_state[0]._Z = MagpieUtils::getZFromZAID(ffzaid);
  initial_state[0]._mass = MagpieUtils::getAFromZAID(ffzaid);

  // sample neutrons per fission
  auto neutrons_per_fission = sampleNu(energy, key);

  // mass balance to find Z and A of second fission product
  initial_state[1]._Z = target_Z - initial_state[0]._Z;
  initial_state[1]._mass = target_A - initial_state[0]._mass - neutrons_per_fission;

  // calculate the kinetic energy of each fission product
  initial_state[0]._energy = total_ke / (1.0 + (initial_state[0]._mass / initial_state[1]._mass));
  initial_state[1]._energy = (initial_state[0]._mass / initial_state[1]._mass) * initial_state[0]._energy;

  // uniformly sample mu and phi (mu is backwards I think but it doesn't matter for uniform sampling)
  // angles[0] = mu, angles[1] = phi
  auto angles = sampleUniformDirection();

  // calculate direction of fission products
  // x: direction(0), y: direction(1), z: direction(2)
  initial_state[0]._direction(0) = std::sqrt(1.0 - angles[0] * angles[0]) * std::cos(angles[1]);
  initial_state[0]._direction(1) = std::sqrt(1.0 - angles[0] * angles[0]) * std::sin(angles[1]);
  initial_state[0]._direction(2) = angles[0];
  initial_state[1]._direction(0) = -initial_state[0]._direction(0);
  initial_state[1]._direction(1) = -initial_state[0]._direction(1);
  initial_state[1]._direction(2) = -initial_state[0]._direction(2);
}

void
DiscreteFissionPKAPDF::readFissionData(const std::vector<unsigned int> & ZAID)
{
  _fission_zaids.resize(4);
  _fission_cdf.resize(4);
  MagpieUtils::neutronEnergyTypes etypes;
  for (unsigned int j = MagpieUtils::Thermal; j != MagpieUtils::enum_max; ++j)
  {
    //Dummy maps
    std::map<unsigned int, std::vector<unsigned int> > zaid_map;
    std::map<unsigned int, std::vector<Real> > cdf_map;
    for (unsigned int i = 0; i < ZAID.size(); ++i)
    {
      std::string the_current_enum = MagpieUtils::neutronEnergyString[j];
      auto path = std::getenv("ENDF_FP_DIR");

      //check if $ENDF_FP_DIR is set
      if (path == NULL)
        mooseError("Set $ENDF_FP_DIR to the directory holding the sum yield data files.");

      std::string ss = path + std::to_string(ZAID[i]) + "_" + the_current_enum + ".txt";

      // check if file exists, if not continue
      auto x = ss;
      const char * filename = x.c_str();
      bool result = MooseUtils::checkFileReadable(filename, false, false);
      if (result == false)
         continue;

      // read the data
      std::ifstream infile(filename);
      std::vector<unsigned int> _zaid_target;
      std::vector<Real> _fission_probabilities;
      int aa;
      Real bb;
      Real prob_prev = 0.0;
      //Accumulates CDF as it reads from the file
      while (infile >> aa >> bb)
      {
        _zaid_target.push_back(aa);
        _fission_probabilities.push_back(bb + prob_prev);
        prob_prev += bb;
      }

      //Renormalize to ensure that cdf[-1] == 1
      unsigned int last_index = _fission_probabilities.size() - 1;
      Real last_value = _fission_probabilities[last_index];
      for (auto & zaid: _fission_probabilities)
         zaid /= last_value;

      //Step 4: Store std::vectors into maps
      zaid_map[ZAID[i]] = _zaid_target;
      cdf_map[ZAID[i]] = _fission_probabilities;
    }
    //Store map into vector for given energy type
    _fission_zaids[j] = zaid_map;
    _fission_cdf[j] = cdf_map;
  }
}

Real
DiscreteFissionPKAPDF::determineFragmentsEnergy(unsigned int Z, unsigned int A)
{
  return 0.1178 * (std::pow(Z, 2.0) / std::pow(A, 1.0 / 3.0)) + 5.8;
}

/**
 * Given the target isotope and neutron energy,
 * determine the avg. number of neutrons per fission.
 */
unsigned int
DiscreteFissionPKAPDF::sampleNu(MagpieUtils::neutronEnergyTypes energy_type, unsigned int zaid)
{
  Real nu_bar;
  if (zaid == 922350 && energy_type == MagpieUtils::Thermal)
    nu_bar = 2.4355;
  else if (zaid == 922330 && energy_type == MagpieUtils::Epithermal)
    nu_bar = 2.4968;
  else if (zaid == 902320 && energy_type == MagpieUtils::Epithermal)
    nu_bar = 2.456;
  else if (zaid == 942390 && energy_type == MagpieUtils::Thermal)
    nu_bar = 2.8836;
  else if (zaid == 942390 && energy_type == MagpieUtils::Epithermal)
    nu_bar = 2.8836;
  else if (zaid == 942400 && energy_type == MagpieUtils::Fast)
    nu_bar = 3.086;
  else if (zaid == 942410 && energy_type == MagpieUtils::Thermal)
    nu_bar = 2.4979;
  else if (zaid == 942420 && energy_type == MagpieUtils::Epithermal)
    nu_bar = 3.189;
  else if (zaid == 952410 && energy_type == MagpieUtils::Thermal)
    nu_bar = 3.239;
  else if (zaid == 962420 && energy_type == MagpieUtils::Epithermal)
    nu_bar = 2.529;
  else if (zaid == 962430 && energy_type == MagpieUtils::Thermal)
    nu_bar = 3.433;
  else if (zaid == 962430 && energy_type == MagpieUtils::Epithermal)
    nu_bar = 3.433;
  else
    nu_bar = 2.5;

  if (nu_bar <= std::floor(nu_bar) + MooseRandom::rand())
    return std::floor(nu_bar);
  else
    return std::ceil(nu_bar);
}
