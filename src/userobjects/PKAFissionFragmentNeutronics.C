/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAFissionFragmentNeutronics.h"
#include "PKAGeneratorBase.h"
#include "PKAGeneratorNeutronicsBase.h"
#include "MultiIndex.h"

registerMooseObject("MagpieApp", PKAFissionFragmentNeutronics);

template <>
InputParameters
validParams<PKAFissionFragmentNeutronics>()
{
  InputParameters params = validParams<PKAGeneratorNeutronicsBase>();
  params.addClassDescription(
      "PKA generator (fission) user object.\n Takes pdf and samples PKAs due to fission.");
  return params;
}

PKAFissionFragmentNeutronics::PKAFissionFragmentNeutronics(const InputParameters & parameters)
  : PKAGeneratorNeutronicsBase(parameters)
{
}

void
PKAFissionFragmentNeutronics::setPDF(const std::vector<unsigned int> & ZAID,
                                     const std::vector<Real> & energies,
                                     const MultiIndex<Real> & probabilities)
{
  _pdf = DiscreteFissionPKAPDF(ZAID, energies, probabilities);
}

void
PKAFissionFragmentNeutronics::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                                         const MyTRIMRasterizer::PKAParameters & pka_parameters,
                                         const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  const auto dt = pka_parameters._dt;
  const auto vol = pka_parameters._volume;
  const auto recoil_rate_scaling = pka_parameters._recoil_rate_scaling;

  mooseAssert(dt >= 0,
              "Passed a negative time window into PKAFissionFragmentNeutronics::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAFissionFragmentNeutronics::appendPKAs");

  if (averaged_data._elements.size() != _partial_neutronics_reaction_rates.size())
    mooseError("Size of averaged_data and partial_reaction_rates must be equal");

  for (unsigned int nuclide = 0; nuclide < _partial_neutronics_reaction_rates.size(); ++nuclide)
  {
    unsigned int num_fission =
        std::floor(recoil_rate_scaling * dt * vol * (*_partial_neutronics_reaction_rates[nuclide]) /
                       (*_averaged_number_densities[nuclide] * averaged_data._site_volume) *
                       averaged_data._elements[nuclide] +
                   getRandomReal());

    for (unsigned i = 0; i < num_fission; ++i)
    {
      std::vector<MyTRIM_NS::IonBase> ion;
      // at this point sample will have Z, m, E
      _pdf.drawSample(ion);

      // set stopping criteria
      ion[0].setEf();
      ion[1].setEf();

      // the tag is the element this PKA get registered as upon stopping
      // -1 means the PKA will be ignored
      ion[0]._tag = ionTag(pka_parameters, ion[0]._Z, ion[0]._m);
      ion[1]._tag = ionTag(pka_parameters, ion[1]._Z, ion[1]._m);

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
