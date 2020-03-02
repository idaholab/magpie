/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MDGranularPropertyAux.h"
#include "MDRunBase.h"

registerMooseObject("MagpieApp", MDGranularPropertyAux);

InputParameters
MDGranularPropertyAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("user_object", "Name of MD runner UserObject");
  params.addParam<MultiMooseEnum>("md_particle_property",
                                  MDRunBase::mdParticleProperties(),
                                  "Property that is injected into auxiliary variable.");
  params.addParam<MooseEnum>("average_type",
                             MDGranularPropertyAux::mdAveragingType(),
                             "The type of average to be taken: "
                             "granular_sum|granular_densitygranular_interstitial_density.");
  params.addClassDescription(
      "Injects properties collected for MD particles from MDRunBase object user_object "
      "into auxiliary variable.");
  return params;
}

MDGranularPropertyAux::MDGranularPropertyAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _md_uo(getUserObject<MDRunBase>("user_object")),
    _average_type(getParam<MooseEnum>("average_type"))
{
  // check length of MultiMooseEnum parameter, get it and check that UO has it
  MultiMooseEnum mme = getParam<MultiMooseEnum>("md_particle_property");
  if (mme.size() != 1)
    mooseError("md_particle_property must contain a single property.");
  _property_id = mme.get(0);
  if (!_md_uo.properties().contains(mme))
    mooseError("Property ", _property_id, " not available from user_object.");

  // ensure MD particles are granular
  if (!_md_uo.isGranular())
    mooseError("user_object stores non-granular particles.");

  // ensure variable is elemental
  if (isNodal())
    mooseError("MDGranularPropertyAux only permits elemental variables.");
}

Real
MDGranularPropertyAux::computeValue()
{
  if (_qp == 0)
  {
    _property_value = 0.0;

    // get the overlapping MD particles
    std::vector<std::pair<unsigned int, Real>> gran_vol;
    _md_uo.granularElementVolumes(_current_elem->unique_id(), gran_vol);

    // loop over the overlapping MD particles and add property value
    Real denominator = 0;
    for (auto & p : gran_vol)
    {
      _property_value += p.second * _md_uo.particleProperty(p.first, _property_id);
      denominator += p.second;
    }

    // compute property value depending on what average type is requested
    if (_average_type == 1)
      _property_value /= _current_elem->volume();
    else if (_average_type == 2)
    {
      if (denominator == 0.0)
        _property_value = 0.0;
      else
        _property_value /= denominator;
    }
  }
  return _property_value;
}

MooseEnum
MDGranularPropertyAux::mdAveragingType()
{
  return MooseEnum("granular_sum=0 granular_density=1 granular_interstitial_density=2");
}
