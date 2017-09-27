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

template<>
InputParameters validParams<PKAFissionFragmentNeutronics>()
{
  InputParameters params = validParams<PKAGeneratorNeutronicsBase>();
  params.addParam<std::vector<PostprocessorName>>("partial_fission_rates", "Partial fission rates per unit volume [sum_g sigma_{f,g,i} * phi_g]. "
                                                  "Provide number density as variable in rasterizer!");
  params.addClassDescription("PKA generator (fission) user object.\n Takes pdf and samples PKAs due to fission.");
  return params;
}

PKAFissionFragmentNeutronics::PKAFissionFragmentNeutronics(const InputParameters & parameters):
    PKAGeneratorNeutronicsBase(parameters)
{
  if (isParamValid("partial_fission_rates"))
  {
    std::vector<PostprocessorName> names = getParam<std::vector<PostprocessorName>>("partial_fission_rates");
    _partial_fission_rates.resize(names.size());
    _stored_pps.resize(names.size());
    for (unsigned int j = 0; j < names.size(); ++j)
      if (_fe_problem.hasPostprocessor(names[j]))
        _partial_fission_rates[j] = &getPostprocessorValueByName(names[j]);
      else
      {
        Real real_value = -std::numeric_limits<Real>::max();
        std::istringstream ss(names[j]);

        if (ss >> real_value && ss.eof())
          _stored_pps[j] = real_value;
        else
          mooseError("Illegal entry in partial_fission_rates: ", names[j]);

        _partial_fission_rates[j] = &_stored_pps[j];
      }
  }
  else
  {
    _stored_pps = {1.0e-8};
    _partial_fission_rates = {&_stored_pps[0]};
  }
}

PKAFissionFragmentNeutronics::~PKAFissionFragmentNeutronics()
{
}

void
PKAFissionFragmentNeutronics::setPDF(const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities)
{
  _pdf = DiscreteFissionPKAPDF(ZAID, energies, probabilities);
}

void
PKAFissionFragmentNeutronics::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  mooseAssert(dt >= 0, "Passed a negative time window into PKAFissionFragmentNeutronics::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAFissionFragmentNeutronics::appendPKAs");

  if (averaged_data._elements.size() != _partial_fission_rates.size())
    mooseError("Size of averaged_data and partial_fission_rates must be equal");

  for (unsigned int nuclide = 0; nuclide < _partial_fission_rates.size(); ++nuclide)
  {
    unsigned int num_fission = std::floor(dt * vol * (*_partial_fission_rates[nuclide]) * averaged_data._elements[nuclide] + getRandomReal());

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
      ion[0]._tag = -1;
      ion[1]._tag = -1;

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
