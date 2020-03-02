/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAGeneratorRecoil.h"
#include "DiscreteFissionPKAPDF.h"
#include "MultiIndex.h"

registerMooseObject("MagpieApp", PKAGeneratorRecoil);

InputParameters
PKAGeneratorRecoil::validParams()
{
  InputParameters params = PKAGeneratorNeutronicsBase::validParams();
  params.addClassDescription(
      "PKA recoil generator user object.\n Takes pdf and samples PKAs due to recoil reaction.");
  return params;
}

PKAGeneratorRecoil::PKAGeneratorRecoil(const InputParameters & parameters)
  : PKAGeneratorNeutronicsBase(parameters)
{
}

void
PKAGeneratorRecoil::setPDF(const std::vector<unsigned int> & ZAID,
                           const std::vector<Real> & energies,
                           const MultiIndex<Real> & probabilities)
{
  _pdf = DiscretePKAPDF(ZAID, energies, probabilities);
#if DEBUG
  _console << "\nPKAGeneratorRecoil object received the following pdf from DiscretePKAPDF object:\n"
           << _pdf;
#endif
}

void
PKAGeneratorRecoil::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                               const MyTRIMRasterizer::PKAParameters & pka_parameters,
                               const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  const auto dt = pka_parameters._dt;
  const auto vol = pka_parameters._volume;
  const auto recoil_rate_scaling = pka_parameters._recoil_rate_scaling;

  mooseAssert(dt >= 0, "Passed a negative time window into PKAGeneratorRecoil::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAGeneratorRecoil::appendPKAs");

  if (averaged_data._elements.size() != _partial_neutronics_reaction_rates.size())
    mooseError("Size of averaged_data and partial_reaction_rates must be equal");

  for (unsigned int nuclide = 0; nuclide < _partial_neutronics_reaction_rates.size(); ++nuclide)
  {
    unsigned int num_recoils =
        std::floor(recoil_rate_scaling * dt * vol * (*_partial_neutronics_reaction_rates[nuclide]) /
                       (*_averaged_number_densities[nuclide]) * averaged_data._elements[nuclide] +
                   getRandomReal());

    for (unsigned i = 0; i < num_recoils; ++i)
    {
      std::vector<MyTRIM_NS::IonBase> ion;
      _pdf.drawSample(ion);
      ion[0].setEf();
      // we need to track this recoil for getting interstitial count right
      ion[0]._tag = nuclide;
      setPosition(ion[0]);
      ion_list.push_back(ion[0]);
    }
  }
}
