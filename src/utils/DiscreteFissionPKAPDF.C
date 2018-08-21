/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "DiscreteFissionPKAPDF.h"
#include "MooseError.h"
#include "MooseRandom.h"
#include "MagpieUtils.h"
#include "MooseUtils.h"

#include "libmesh/utility.h"

#include <fstream>
#include <map>

DiscreteFissionPKAPDF::DiscreteFissionPKAPDF()
  : DiscretePKAPDFBase(),
    _marginal_cdf_target(MultiIndex<Real>({1})),
    _conditional_cdf_energy(MultiIndex<Real>({1}))
{
}

DiscreteFissionPKAPDF::DiscreteFissionPKAPDF(const std::vector<unsigned int> & ZAID,
                                             const std::vector<Real> & energies,
                                             const MultiIndex<Real> & probabilities)
  : DiscretePKAPDFBase(ZAID, energies),
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
  computeMagnitude(probabilities);
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
  for (auto target : _marginal_cdf_target)
  {
    index = target.first;
    index[0] -= 1;
    if (target.first[0] > 0)
      target.second += _marginal_cdf_target(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  index[0] = _nZA - 1;
  Real last_value = _marginal_cdf_target(index);
  for (auto target : _marginal_cdf_target)
    target.second /= last_value;

  // Step 1: Compute the conditional distribution of the target w.r.t the sampled energy
  _conditional_cdf_energy = probabilities;

  // Step 2: Compute cdf
  for (auto energy : _conditional_cdf_energy)
  {
    index = energy.first;
    index[1] -= 1;
    if (energy.first[1] > 0)
      energy.second += _conditional_cdf_energy(index);
  }

  // Step 3: Renormalize to ensure that cdf[-1] == 1
  std::vector<Real> last_element_per_target;
  for (unsigned int j = 0; j < _conditional_cdf_energy.size()[0]; ++j)
    last_element_per_target.push_back(_conditional_cdf_energy({j, _ng - 1}));
  for (auto energy : _conditional_cdf_energy)
    energy.second /= last_element_per_target[energy.first[0]];
}

void
DiscreteFissionPKAPDF::drawSample(std::vector<MyTRIM_NS::IonBase> & initial_state) const
{
  // resize initial_state
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
  auto neutron_energy =
      (_energies[sampled_indices[1]] - _energies[sampled_indices[1] + 1]) * MooseRandom::rand() +
      _energies[sampled_indices[1] + 1];

  // pull the right fission yield table and sample Z, A, and energy
  auto energy = MagpieUtils::determineNeutronType(neutron_energy);
  auto key = _zaids[sampled_indices[0]];

  // pull the map based on the neutron energy
  auto zaid_map = _fission_zaids[energy];
  auto cdf_map = _fission_cdf[energy];

  /// check if key exists
  if (zaid_map.find(key) == zaid_map.end())
    mooseError("Fission product ZAID's not found for this target isotope, ",
               key,
               " at this energy, ",
               energy);

  if (cdf_map.find(key) == cdf_map.end())
    mooseError("Sum yield not found for this target isotope, ", key, " at this energy, ", energy);

  // now we know the key1/key2 pair exists
  auto zaid_products = zaid_map[key];
  auto cdf_products = cdf_map[key];

  // calculate the total kinetic energy of the fission products
  auto total_ke = determineFragmentsEnergy(target_Z, target_A);

  // sample the fission yield cdf to determine the first fission product
  auto ffindex = sampleHelper(cdf_products);
  auto ffzaid = zaid_products[ffindex];
  initial_state[0]._Z = MagpieUtils::getZFromZAID(ffzaid);
  initial_state[0]._m = MagpieUtils::getAFromZAID(ffzaid);

  // sample neutrons per fission
  auto neutrons_per_fission = sampleNu(energy, key);

  // mass balance to find Z and A of second fission product
  initial_state[1]._Z = target_Z - initial_state[0]._Z;
  initial_state[1]._m = target_A - initial_state[0]._m - neutrons_per_fission;

  // calculate the kinetic energy of each fission product
  initial_state[0]._E = total_ke / (1.0 + (initial_state[0]._m / initial_state[1]._m));
  initial_state[1]._E = (initial_state[0]._m / initial_state[1]._m) * initial_state[0]._E;
}

void
DiscreteFissionPKAPDF::readFissionData(const std::vector<unsigned int> & zaid_list)
{
  _fission_zaids.resize(MagpieUtils::NET_MAX);
  _fission_cdf.resize(MagpieUtils::NET_MAX);

  // iterate over all neutron energy types
  for (unsigned int energy = 0; energy < MagpieUtils::NET_MAX; ++energy)
  {
    // Dummy maps
    std::map<unsigned int, std::vector<unsigned int>> zaid_map;
    std::map<unsigned int, std::vector<Real>> cdf_map;
    for (auto & zaid : zaid_list)
    {
      // check if $ENDF_FP_DIR is set
      auto path = std::getenv("ENDF_FP_DIR");
      if (path == NULL)
        mooseError("Set $ENDF_FP_DIR to the directory holding the sum yield data files.");

      // check if file exists, if not continue
      std::string filename =
          path + std::to_string(zaid) + "_" + MagpieUtils::neutronEnergyName(energy) + ".txt";
      if (!MooseUtils::checkFileReadable(filename, false, false))
      {
        // most isotopes have fission data for High but not for Fast wich is weird
        // if fast is not available, then try to load High instead
        if (MagpieUtils::neutronEnergyName(energy) == "Fast")
        {
          std::string alt_filename = path + std::to_string(zaid) + "_" + "High" + ".txt";
          if (!MooseUtils::checkFileReadable(alt_filename, false, false))
            continue;
          filename = alt_filename;
        }
        else
          continue;
      }

      // read the data
      std::ifstream infile(filename.c_str());
      std::vector<unsigned int> zaid_target;
      std::vector<Real> fission_probabilities;
      int aa;
      Real bb;
      Real prob_prev = 0.0;
      while (infile >> aa >> bb)
      {
        zaid_target.push_back(aa);
        fission_probabilities.push_back(bb + prob_prev);

        // accumulate CDF as it reads from the file
        prob_prev += bb;
      }

      // renormalize to ensure that cdf[-1] == 1
      unsigned int last_index = fission_probabilities.size() - 1;
      Real last_value = fission_probabilities[last_index];
      for (auto & probability : fission_probabilities)
        probability /= last_value;

      // Step 4: Store std::vectors into maps
      zaid_map[zaid] = zaid_target;
      cdf_map[zaid] = fission_probabilities;
    }

    // store map into vector for given energy type
    _fission_zaids[energy] = zaid_map;
    _fission_cdf[energy] = cdf_map;
  }
}

Real
DiscreteFissionPKAPDF::determineFragmentsEnergy(unsigned int Z, unsigned int A) const
{
  // MyTRIM expects energies to be in eV, hence the 1.0e6 factor
  return (0.1178 * (Utility::pow<2>(Z) / std::cbrt(A)) + 5.8) * 1.0e6;
}

/**
 * Given the target isotope and neutron energy,
 * determine the avg. number of neutrons per fission.
 */
unsigned int
DiscreteFissionPKAPDF::sampleNu(MagpieUtils::NeutronEnergyType energy_type, unsigned int zaid) const
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

void
DiscreteFissionPKAPDF::computeMagnitude(MultiIndex<Real> probabilities)
{
  _magnitude = 0.0;
  for (auto it : probabilities)
  {
    MultiIndex<Real>::size_type index = it.first;
    Real delE = _energies[index[1]] - _energies[index[1] + 1];
    _magnitude += delE * it.second;
  }
}
