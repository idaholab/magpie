/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAGeneratorNeutronicsBase.h"
#include "MagpieUtils.h"
#include "MultiIndex.h"
#include "DiscreteFissionPKAPDF.h"

template<>
InputParameters validParams<PKAGeneratorNeutronicsBase>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addParam<std::vector<PostprocessorName>>("partial_reaction_rates", "Partial neutronic reaction rates per unit volume [sum_g xs_{r,g,i} * phi_g], "
                                                  "r: reaction type, g: energy group, i: nuclide id.Provide number density as variable in rasterizer!");
  params.addParam<std::vector<PostprocessorName>>("averaged_number_densities", "The number density of the species averaged over the domain.");
  params.addClassDescription("PKA generator (neutronics) user object base class.\n Takes pdf and samples PKAs due to various interactions.");
  return params;
}

PKAGeneratorNeutronicsBase::PKAGeneratorNeutronicsBase(const InputParameters & parameters) :
    PKAGeneratorBase(parameters)
{
  if (isParamValid("partial_reaction_rates"))
  {
    std::vector<PostprocessorName> names = getParam<std::vector<PostprocessorName>>("partial_reaction_rates");
    _partial_neutronics_reaction_rates.resize(names.size());
    _stored_reaction_rates.resize(names.size());
    for (unsigned int j = 0; j < names.size(); ++j)
      if (_fe_problem.hasPostprocessor(names[j]))
        _partial_neutronics_reaction_rates[j] = &getPostprocessorValueByName(names[j]);
      else
      {
        Real real_value = -std::numeric_limits<Real>::max();
        std::istringstream ss(names[j]);

        if (ss >> real_value && ss.eof())
          _stored_reaction_rates[j] = real_value;
        else
          mooseError("Illegal entry in _partial_neutronics_reaction_rates: ", names[j]);

        _partial_neutronics_reaction_rates[j] = &_stored_reaction_rates[j];
      }
  }
  else
  {
    _stored_reaction_rates = {1.0e-8};
    _partial_neutronics_reaction_rates = {&_stored_reaction_rates[0]};
  }

  if (isParamValid("averaged_number_densities"))
  {
    std::vector<PostprocessorName> names = getParam<std::vector<PostprocessorName>>("averaged_number_densities");
    _averaged_number_densities.resize(names.size());
    _stored_densities.resize(names.size());
    for (unsigned int j = 0; j < names.size(); ++j)
      if (_fe_problem.hasPostprocessor(names[j]))
        _averaged_number_densities[j] = &getPostprocessorValueByName(names[j]);
      else
      {
        Real real_value = -std::numeric_limits<Real>::max();
        std::istringstream ss(names[j]);

        if (ss >> real_value && ss.eof())
          _stored_densities[j] = real_value;
        else
          mooseError("Illegal entry in _partial_neutronics_reaction_rates: ", names[j]);

        _averaged_number_densities[j] = &_stored_reaction_rates[j];
      }
  }
  else
  {
    _averaged_number_densities.resize(_stored_reaction_rates.size());
    _stored_densities.resize(_stored_reaction_rates.size());
    for (unsigned int j = 0; j < _stored_densities.size(); ++j)
    {
      _stored_densities[j] = 1;
      _averaged_number_densities[j] = &_stored_densities[j];
    }
  }

  if (_averaged_number_densities.size() != _partial_neutronics_reaction_rates.size())
    mooseError("partial_reaction_rates and averaged_number_densities must have the same number of entries.");
}
