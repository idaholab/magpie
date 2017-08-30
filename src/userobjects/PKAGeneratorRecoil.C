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

template<>
InputParameters validParams<PKAGeneratorRecoil>()
{
  InputParameters params = validParams<PKAGeneratorNeutronicsBase>();
  params.addParam<std::vector<PostprocessorName>>("partial_recoil_rates", "Partial recoil rates per unit volume [recoil reaction rate per nuclide]. "
                                                  "Provide number density as variable in rasterizer!");
  params.addClassDescription("PKA recoil generator user object.\n Takes pdf and samples PKAs due to recoil reaction.");
  return params;
}

PKAGeneratorRecoil::PKAGeneratorRecoil(const InputParameters & parameters):
    PKAGeneratorNeutronicsBase(parameters)
{
  if (isParamValid("partial_recoil_rates"))
  {
    std::vector<PostprocessorName> names = getParam<std::vector<PostprocessorName>>("partial_recoil_rates");
    _partial_recoil_rates.resize(names.size());
    _stored_pps.resize(names.size());
    for (unsigned int j = 0; j < names.size(); ++j)
      if (_fe_problem.hasPostprocessor(names[j]))
        _partial_recoil_rates[j] = &getPostprocessorValueByName(names[j]);
      else
      {
        Real real_value = -std::numeric_limits<Real>::max();
        std::istringstream ss(names[j]);

        if (ss >> real_value && ss.eof())
          _stored_pps[j] = real_value;
        else
          mooseError("Illegal entry in partial_recoil_rates: ", names[j]);

        _partial_recoil_rates[j] = &_stored_pps[j];
      }
  }
  else
  {
    _stored_pps = {1.0e-8};
    _partial_recoil_rates = {&_stored_pps[0]};
  }
}

void
PKAGeneratorRecoil::setPDF(const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities)
{
  _pdf = DiscretePKAPDF(ZAID, energies, probabilities);
}

void
PKAGeneratorRecoil::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  mooseAssert(dt >= 0, "Passed a negative time window into PKAGeneratorRecoil::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAGeneratorRecoil::appendPKAs");

  if (averaged_data._elements.size() != _partial_recoil_rates.size())
    mooseError("Size of averaged_data and partial_recoil_rates must be equal");

  for (unsigned int nuclide = 0; nuclide < _partial_recoil_rates.size(); ++nuclide)
  {
    unsigned int num_recoils = std::floor(dt * vol * (*_partial_recoil_rates[nuclide]) * averaged_data._elements[nuclide] + getRandomReal());

    for (unsigned i = 0; i < num_recoils; ++i)
    {
      std::vector<MyTRIM_NS::IonBase> ion;
      _pdf.drawSample(ion);
      ion[0].setEf();
      // the tag is the element this PKA get registered as upon stopping
      // -1 means the PKA will be ignored
      ion[0]._tag = -1;
      // set location of the fission event
      setPosition(ion[0]);
      // add PKAs to list
      ion_list.push_back(ion[0]);
    }
  }
}
