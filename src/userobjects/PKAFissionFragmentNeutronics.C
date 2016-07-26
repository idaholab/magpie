#include "PKAFissionFragmentNeutronics.h"
#include "PKAGeneratorBase.h"
#include "PKAGeneratorNeutronicsBase.h"
#include "MultiIndex.h"

template<>
InputParameters validParams<PKAFissionFragmentNeutronics>()
{
  InputParameters params = validParams<PKAGeneratorNeutronicsBase>();
  return params;
}

PKAFissionFragmentNeutronics::PKAFissionFragmentNeutronics(const InputParameters & parameters):
    PKAGeneratorNeutronicsBase(parameters)
{
}

PKAFissionFragmentNeutronics::~PKAFissionFragmentNeutronics()
{
}

void
PKAFissionFragmentNeutronics::setPDF(Real magnitude, const std::vector<unsigned int> & ZAID, const std::vector<Real> & energies, const MultiIndex<Real> & probabilities)
{
  _pdf = DiscreteFissionPKAPDF(magnitude, ZAID, energies, probabilities);
}

void
PKAFissionFragmentNeutronics::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData &) const
{
  mooseAssert(dt > 0, "Passed a negative time window into PKAFissionFragmentNeutronics::appendPKAs");
  mooseAssert(vol > 0, "Passed a volume into PKAFissionFragmentNeutronics::appendPKAs");

  //FIXME need getMagnitude from user object currently in another commit
  //unsigned int num_fission = std::floor(dt * vol * _pdf.getMagnitude() + getRandomReal());
  unsigned int num_fission = 1e5;

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
    ion[0].tag = -1;
    ion[1].tag = -1;

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
