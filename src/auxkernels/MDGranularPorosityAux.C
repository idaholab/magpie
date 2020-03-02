/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MDGranularPorosityAux.h"
#include "MDRunBase.h"

registerMooseObject("MagpieApp", MDGranularPorosityAux);

InputParameters
MDGranularPorosityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("user_object", "Name of MD runner UserObject");
  params.addParam<bool>(
      "compute_packing_fraction", false, "Whether to compute porosity or packing fraction.");
  params.addClassDescription(
      "Computes porosity or packing fraction (1 - porosity) and injects it into and aux variable.");
  return params;
}

MDGranularPorosityAux::MDGranularPorosityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _md_uo(getUserObject<MDRunBase>("user_object")),
    _compute_packing(getParam<bool>("compute_packing_fraction"))
{
  // ensure MD particles are granular
  if (!_md_uo.isGranular())
    mooseError("user_object stores non-granular particles.");

  // ensure variable is elemental
  if (isNodal())
    mooseError("MDGranularPorosityAux only permits elemental variables.");
}

Real
MDGranularPorosityAux::computeValue()
{
  if (_qp == 0)
  {
    _packing_fraction = 0.0;

    // get the overlapping MD particles
    std::vector<std::pair<unsigned int, Real>> gran_vol;
    _md_uo.granularElementVolumes(_current_elem->unique_id(), gran_vol);

    // add the overlapping volumes
    for (auto & p : gran_vol)
      _packing_fraction += p.second;

    // divide by element volume
    _packing_fraction /= _current_elem->volume();
  }
  if (_compute_packing)
    return _packing_fraction;
  return 1 - _packing_fraction;
}
