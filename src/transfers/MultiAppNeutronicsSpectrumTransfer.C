#include "MultiAppNeutronicsSpectrumTransfer.h"
#include "MultiApp.h"
#include "FEProblem.h"
#include "MooseTypes.h"
#include "MultiIndex.h"
#include "UserObject.h"
#include "NeutronicsSpectrumSamplerBase.h"
#include "PKAGeneratorNeutronicsBase.h"

template<>
InputParameters validParams<MultiAppNeutronicsSpectrumTransfer>()
{
  InputParameters params = validParams<MultiAppTransfer>();
  params.addRequiredParam<UserObjectName>("pka_neutronics", "PKA generator object name.");
  params.addRequiredParam<UserObjectName>("radiation_damage_sampler", "Neutronics user object providing the PDF data.");
  return params;
}

MultiAppNeutronicsSpectrumTransfer::MultiAppNeutronicsSpectrumTransfer(const InputParameters & parameters) :
    MultiAppTransfer(parameters),
    _pka_generator_name(getParam<UserObjectName>("pka_neutronics")),
    _neutronics_pdf_name(getParam<UserObjectName>("radiation_damage_sampler"))
{
  if (_direction != TO_MULTIAPP)
    mooseError("MultiAppNeutronicsSpectrumTransfer can only send data from a neutronics master app to a mesoscale multiapp.");
}

void
MultiAppNeutronicsSpectrumTransfer::execute()
{
  // get the neutronics PDF user object
  FEProblem & from_problem = _multi_app->problem();
  const NeutronicsSpectrumSamplerBase & neutronics_pdf = from_problem.getUserObject<NeutronicsSpectrumSamplerBase>(_neutronics_pdf_name);

  // loop over all sub apps and copy over the neutronics data
  for (unsigned int i = 0; i < _multi_app->numGlobalApps(); ++i)
    if (_multi_app->hasLocalApp(i))
    {
      std::vector<unsigned int> zaids = neutronics_pdf.getZAIDs();
      std::vector<Real> energies = neutronics_pdf.getEnergies();
      MultiIndex<Real> probabilities = neutronics_pdf.getPDF(i);

      for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
      {
        PKAGeneratorNeutronicsBase & pka_uo = const_cast<PKAGeneratorNeutronicsBase &>(_multi_app->appProblem(i).getUserObject<PKAGeneratorNeutronicsBase>(_pka_generator_name, tid));
        pka_uo.setPDF(zaids, energies, probabilities);
      }
    }
}
