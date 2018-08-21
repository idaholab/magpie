/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAFissionFragmentEmpirical.h"

registerMooseObject("MagpieApp", PKAFissionFragmentEmpirical);

template <>
InputParameters
validParams<PKAFissionFragmentEmpirical>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addClassDescription("Fission fragment PKA generator usimg an empirical mass and energy "
                             "distribution for LWR UO2 fuel");
  params.addParam<PostprocessorName>("fission_rate",
                                     1e-8,
                                     "Fission rate per unit volume (uses mesh units defined in the "
                                     "rasterizer and moose time units)");
  params.addRequiredCoupledVar("relative_density",
                               "Relative UO2 density (1 is fully dense, 0 is no UO2");
  return params;
}

PKAFissionFragmentEmpirical::PKAFissionFragmentEmpirical(const InputParameters & parameters)
  : PKAGeneratorBase(parameters),
    _fission_rate(getPostprocessorValue("fission_rate")),
    _relative_density(coupledValue("relative_density"))
{
}

void
PKAFissionFragmentEmpirical::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                                        const MyTRIMRasterizer::PKAParameters & pka_parameters,
                                        const MyTRIMRasterizer::AveragedData &) const
{
  const auto dt = pka_parameters._dt;
  const auto vol = pka_parameters._volume;
  const auto recoil_rate_scaling = pka_parameters._recoil_rate_scaling;

  mooseAssert(dt >= 0,
              "Passed a negative time window into PKAFissionFragmentEmpirical::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAFissionFragmentEmpirical::appendPKAs");

  unsigned int num_fission = std::floor(
      recoil_rate_scaling * dt * vol * _fission_rate * _relative_density[0] + getRandomReal());

  /// Mass inverter to sample PKA mass distribution
  MyTRIM_NS::MassInverter mass_inverter;

  for (unsigned i = 0; i < num_fission; ++i)
  {
    // each fission event generates a pair of recoils
    MyTRIM_NS::IonBase ion1, ion2;

    // sample fission fragment masses
    ion1._m = mass_inverter.x(getRandomReal());
    ion2._m = 235.0 - ion1._m - 2.0; // thermal fission emits 2n

    // Total energy in eV (energy inverter output MeV)
    MyTRIM_NS::EnergyInverter energy_inverter;
    energy_inverter.setMass(ion1._m);
    Real Etot = energy_inverter.x(getRandomReal()) * 1e6;
    ion1._E = Etot * ion2._m / (ion1._m + ion2._m);
    ion2._E = Etot - ion1._E;

    // assume p/n ratio like U (Semi-empirical mass formula would be marginally better ~15%,
    // but it is harder to achieve charge neutrality
    ion1._Z = std::round((ion1._m * 92.0) / (235.0 - 2.0));
    ion2._Z = 92 - ion1._Z;

    // set stopping criteria
    ion1.setEf();
    ion2.setEf();

    // the tag is the element this PKA get registered as upon stopping
    // -1 means the PKA will be ignored
    ion1._tag = ionTag(pka_parameters, ion1._Z, ion1._m);
    ion2._tag = ionTag(pka_parameters, ion2._Z, ion2._m);

    // set location of the fission event
    setPosition(ion1);
    ion2._pos = ion1._pos;

    // set random direction for ion 1 and opposite direction for ion 2
    setRandomDirection(ion1);
    ion2._dir = -ion1._dir;

    // add PKAs to list
    ion_list.push_back(ion1);
    ion_list.push_back(ion2);
  }
}
