#include "MultiAppRadiationDamageTransfer.h"
#include "MultiApp.h"
#include "FEProblem.h"
#include "MooseTypes.h"
#include "MultiIndex.h"
#include "UserObject.h"
#include "RadiationDamageBase.h"
#include "PKAGeneratorNeutronicsBase.h"

template<>
InputParameters validParams<MultiAppRadiationDamageTransfer>()
{
  InputParameters params = validParams<MultiAppTransfer>();
  params.addRequiredParam<UserObjectName>("pka_neutronics", "PKA generator object name.");
  params.addRequiredParam<UserObjectName>("radiation_damage_sampler", "Neutronics user object providing the PDF data.");
  return params;
}

MultiAppRadiationDamageTransfer::MultiAppRadiationDamageTransfer(const InputParameters & parameters) :
    MultiAppTransfer(parameters),
    _pka_generator_name(getParam<UserObjectName>("pka_neutronics")),
    _neutronics_pdf_name(getParam<UserObjectName>("radiation_damage_sampler"))
{
  if (_direction != TO_MULTIAPP)
    mooseError("MultiAppRadiationDamageTransfer can only send data from a neutronics master app to a mesoscale multiapp.");
}

void
MultiAppRadiationDamageTransfer::execute()
{
  // get the neutronics PDF user object
  FEProblem & from_problem = _multi_app->problem();
  const RadiationDamageBase & neutronics_pdf = from_problem.getUserObject<RadiationDamageBase>(_neutronics_pdf_name);

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
