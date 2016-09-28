#include "PKAFissionFragmentNeutronics.h"
#include "PKAGeneratorBase.h"
#include "PKAGeneratorNeutronicsBase.h"
#include "MultiIndex.h"

template<>
InputParameters validParams<PKAFissionFragmentNeutronics>()
{
  InputParameters params = validParams<PKAGeneratorNeutronicsBase>();
  params.addParam<Real>("fission_rate", 1e-8, "Fission rate per unit volume (default is in nm^-3*s^-1)");
  params.addClassDescription("PKA generator (fission) user object.\n Takes pdf and samples PKAs due to fission.");
  return params;
}

PKAFissionFragmentNeutronics::PKAFissionFragmentNeutronics(const InputParameters & parameters):
    PKAGeneratorNeutronicsBase(parameters),
    _fission_rate(getParam<Real>("fission_rate"))
{
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

  //FIXME this only works for a single source of PKAs dependent on a single variable
  // This variable must be provided as the first var argument in the rasterizer
  if (averaged_data._elements.size() > 1)
    mooseDoOnce(mooseWarning("PKAFissionFragmentNeutronics::appendPKAs only works for a single PDF associated with a single variable in the rasterizer since there is only one fissionable material for now."));

  unsigned int num_fission = std::floor(dt * vol * _fission_rate * averaged_data._elements[0] + getRandomReal());

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
