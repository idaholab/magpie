/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAGeneratorAlphaDecay.h"
#include "MagpieUtils.h"
#include <algorithm>

registerMooseObject("MagpieApp", PKAGeneratorAlphaDecay);

template <>
InputParameters
validParams<PKAGeneratorAlphaDecay>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addClassDescription("A PKAGenerator for starting alpha particles from decay\nDecay data "
                             "is retrieved from file data/alpha_decay/alpha_decay.txt");
  params.addRequiredParam<std::vector<unsigned int>>("ZAID",
                                                     "The Z/A ids corresponding to the var "
                                                     "arguments in the rasterizer [1e4 * A + 10 * "
                                                     "Z + state]");
  MooseEnum timeUnit("second=0 millisecond=1 microsecond=2", "second");
  params.addParam<MooseEnum>("time_unit", timeUnit, "Unit of time used in this model");
  return params;
}

PKAGeneratorAlphaDecay::PKAGeneratorAlphaDecay(const InputParameters & parameters)
  : PKAGeneratorBase(parameters), _zaids(getParam<std::vector<unsigned int>>("ZAID"))
{
  switch (getParam<MooseEnum>("time_unit"))
  {
    case 0:
      _time_conversion = 1.0;
      break;
    case 1:
      _time_conversion = 1.0e3;
      break;
    case 2:
      _time_conversion = 1.0e6;
      break;
  }
  readAlphaData();
}

void
PKAGeneratorAlphaDecay::readAlphaData()
{
  auto path = std::getenv("ALPHA_DIR");
  if (path == NULL)
    mooseError("Set $ALPHA_DIR to the directory holding the alpha_decay.txt file.");

  std::string s = "alpha_decay.txt";
  std::string filename = path + s;

  if (!MooseUtils::checkFileReadable(filename, false, false))
    mooseError("Missing/Corrupted data/alpha_decay/alpha_decay.txt file");

  std::ifstream infile(filename.c_str());
  unsigned int nentries;
  infile >> nentries;
  for (unsigned int j = 0; j < nentries; ++j)
  {
    unsigned int zaid;
    unsigned int n_modes;
    infile >> zaid >> n_modes;
    std::vector<DecayData> decay(n_modes, DecayData());
    for (unsigned l = 0; l < n_modes; ++l)
    {
      unsigned int Z, A;
      Real half_life, energy, intensity;
      infile >> Z >> A >> half_life >> energy >> intensity;
      // adjust half_life unit; file units are
      decay[l]._decay_constants = std::log(2.0) / (half_life * _time_conversion);
      decay[l]._alpha_energies = energy;
      decay[l]._intensities = intensity;
    }

    if (std::find(_zaids.begin(), _zaids.end(), zaid) != _zaids.end())
    {
      if (_decay_data_sets.find(zaid) != _decay_data_sets.end())
        mooseError("Identical Z/A id was provided twice");
      _decay_data_sets[zaid] = decay;
    }
  }
  infile.close();

  // complete the data by dealing with absent nuclides ==> stable
  for (auto & z : _zaids)
  {
    if (_decay_data_sets.find(z) == _decay_data_sets.end())
    {
      std::vector<DecayData> decay(1, DecayData());
      decay[0]._decay_constants = 0.0;
      decay[0]._alpha_energies = 0.0;
      decay[0]._intensities = 1.0;
      _decay_data_sets[z] = decay;
    }
  }

// in debug mode we want some info on what is stored in _decay_data_sets
#if DEBUG
  _console << "\nAlpha decay read from file\n";
  for (auto & d : _decay_data_sets)
  {
    _console << "ZAID: " << d.first << " has " << d.second.size() << " decay channels:\n";
    for (unsigned int j = 0; j < d.second.size(); ++j)
      _console << "lambda = " << d.second[j]._decay_constants << " "
               << "energy = " << d.second[j]._alpha_energies << " "
               << "intensity = " << d.second[j]._intensities << "\n";
  }
  _console << std::endl;
#endif
}

void
PKAGeneratorAlphaDecay::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                                   const MyTRIMRasterizer::PKAParameters & pka_parameters,
                                   const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  const auto dt = pka_parameters._dt;
  const auto vol = pka_parameters._volume;
  const auto recoil_rate_scaling = pka_parameters._recoil_rate_scaling;

  mooseAssert(dt >= 0,
              "Passed a negative time window into PKAFissionFragmentNeutronics::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAFissionFragmentNeutronics::appendPKAs");

  if (averaged_data._elements.size() != _zaids.size())
    mooseError("Size of averaged_data and ZAID must be equal");

  for (unsigned int nuclide = 0; nuclide < _zaids.size(); ++nuclide)
  {
    unsigned int zaid = _zaids[nuclide];

    auto it = _decay_data_sets.find(zaid);

    if (it == _decay_data_sets.end())
      mooseError("ZAID ", zaid, " not found in decay data set");

    auto & decay = it->second;
    for (unsigned int l = 0; l < decay.size(); ++l)
    {
      unsigned int num_decay = std::floor(recoil_rate_scaling * vol * decay[l]._intensities *
                                              averaged_data._elements[nuclide] *
                                              (1.0 - std::exp(-decay[l]._decay_constants * dt)) +
                                          getRandomReal());

      for (unsigned i = 0; i < num_decay; ++i)
      {
        std::vector<MyTRIM_NS::IonBase> ion(2);

        // set stopping criteria
        ion[0].setEf();
        ion[1].setEf();

        // get parent Z/A
        unsigned int parent_Z = MagpieUtils::getZFromZAID(zaid);
        unsigned int parent_A = MagpieUtils::getAFromZAID(zaid);

        // set Z/A values; [0] is alpha, [1] is recoil
        ion[0]._Z = 2;
        ion[0]._m = 4.0;

        ion[1]._Z = parent_Z - 2;
        ion[1]._m = parent_A - 4;

        ion[0]._tag = ionTag(pka_parameters, ion[0]._Z, ion[0]._m);
        ion[1]._tag = ionTag(pka_parameters, ion[1]._Z, ion[1]._m);

        // set the energies; use momentum conservation
        ion[0]._E = decay[l]._alpha_energies;
        ion[1]._E = decay[l]._alpha_energies * ion[0]._m / ion[1]._m;

        // set location of the fission event
        setPosition(ion[0]);
        ion[1]._pos = ion[0]._pos;

        // set random direction for ion 1 and opposite direction for ion 2
        setRandomDirection(ion[0]);
        ion[1]._dir = -ion[0]._dir;

        // add PKAs to list
        ion_list.push_back(ion[0]);
        ion_list.push_back(ion[1]);
      }
    }
  }
}
