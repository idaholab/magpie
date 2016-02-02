#include "PKAFissionFragmentEmpirical.h"

template<>
InputParameters validParams<PKAFissionFragmentEmpirical>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addParam<Real>("fission_rate", 1e-8, "Fission rate per unit volume (default is in nm^-3*s^-1)");
  params.addRequiredCoupledVar("relative_density", "Relative UO2 density (1 is fully dense, 0 is no UO2");
  return params;
}

PKAFissionFragmentEmpirical::PKAFissionFragmentEmpirical(const InputParameters & parameters) :
    PKAGeneratorBase(parameters),
    _fission_rate(getParam<Real>("fission_rate")),
    _relative_density(coupledValue("relative_density"))
{
}

void
PKAFissionFragmentEmpirical::appendPKAs(std::vector<MyTRIM_NS::ionBase> & ion_list, Real dt, Real vol) const
{
  mooseAssert(dt > 0, "Passed a negative time window into PKAFissionFragmentEmpirical::appendPKAs");
  mooseAssert(vol > 0, "Passed a volume into PKAFissionFragmentEmpirical::appendPKAs");

  unsigned int num_fission = std::floor(dt * vol * _fission_rate + getRandomReal());

  /// Mass inverter to sample PKA mass distribution
  MyTRIM_NS::massInverter mass_inverter;

  for (unsigned i = 0; i < num_fission; ++i)
  {
    // each fission event generates a pair of recoils
    MyTRIM_NS::ionBase ion1, ion2;

    // sample fission fragment masses
    ion1.m1 = mass_inverter.x(getRandomReal());
    ion2.m1 = 235.0 - ion1.m1 - 2.0; // thermal fission emits 2n

    // Total energy in eV (energy inverter output MeV)
    MyTRIM_NS::energyInverter energy_inverter;
    energy_inverter.setMass(ion1.m1);
    Real Etot = energy_inverter.x(getRandomReal()) * 1e6;
    ion1.e = Etot * ion2.m1 / (ion1.m1 + ion2.m1);
    ion2.e = Etot - ion1.e;

    // assume p/n ratio like U (Semi-empirical mass formula would be marginally better ~15%,
    // but it is harder to achieve charge neutrality
    ion1.z1 = ::round((ion1.m1 * 92.0) / (235.0 - 2.0));
    ion2.z1 = 92 - ion1.z1;

    // set stopping criteria
    ion1.set_ef();
    ion2.set_ef();

    // set location of the fission event
    setPosition(ion1);
    ion2.pos[0] = ion1.pos[0];
    ion2.pos[1] = ion1.pos[1];
    ion2.pos[2] = ion1.pos[2];

    // set random direction for ion 1 and opposite direction for ion 2
    setRandomDirection(ion1);
    ion2.dir[0] = -ion1.dir[0];
    ion2.dir[1] = -ion1.dir[1];
    ion2.dir[2] = -ion1.dir[2];

    // add PKAs to list
    ion_list.push_back(ion1);
    ion_list.push_back(ion2);
  }
}
