#include "PKAFissionFragmentEmpirical.h"
#include "MathUtils.h"

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

  for (unsigned i = 0; i < num_fission; ++i)
  {
    // each fission event generates a pair of recoils
    MyTRIM_NS::ionBase ion1, ion2;

    // sample fission fragment masses
    ion1.m = _mass_inverter.x(getRandomReal());
    ion2.m = 235.0 - ion1.m;
    e->setMass(A1);

    // Total energy in eV (energy inverter output MeV)
    Real Etot = _energy_inverter.x(getRandomReal()) * 1e6;
    ion1.E = Etot * ion2.m / (ion1.m + ion2.m);
    ion2.E = Etot - ion1.E;

    // assume p/n ratio like Xe, TODO: we can probably do better
    ion1.Z = MathUtils::round((A1 * 92.0) / 235.0);
    ion2.Z = 92 - ion1.Z;

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
