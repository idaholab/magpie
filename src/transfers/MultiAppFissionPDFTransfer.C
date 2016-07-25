#include "MultiAppFissionPDFTransfer.h"
#include "MultiApp.h"
#include "FEProblem.h"
#include "MooseTypes.h"
#include "MultiIndex.h"


/**
 * Dummy objects!
 */
#include "UserObject.h"
// Sebastians Neutronics PDF object
class NeutronicsPDF : public UserObject {
public:
  MultiIndex<Real> getPDF(unsigned int sub_app_id) const { return MultiIndex<Real>({1,1}); }
};
// Pedram's PKA Generator base class for neutronics based PKS generators
class PKAGeneratorNeutronicsBase : public UserObject {
public:
  void setPDF(const MultiIndex<Real> & pdf_data) {}
};


template<>
InputParameters validParams<MultiAppFissionPDFTransfer>()
{
  InputParameters params = validParams<MultiAppTransfer>();
  params.addRequiredParam<UserObjectName>("pka_generator", "PKAFissionFragmentNeutronics PKA generator object name.");
  params.addRequiredParam<UserObjectName>("neutronics_pdf", "Neutronics user object providing the PDF data.");
  return params;
}

MultiAppFissionPDFTransfer::MultiAppFissionPDFTransfer(const InputParameters & parameters) :
    MultiAppTransfer(parameters),
    _pka_generator_name(getParam<UserObjectName>("pka_generator")),
    _neutronics_pdf_name(getParam<UserObjectName>("neutronics_pdf"))
{
  if (_direction != TO_MULTIAPP)
    mooseError("MultiAppFissionPDFTransfer can only send data from a neutronics master app to a mesoscale multiapp.");
}

void
MultiAppFissionPDFTransfer::execute()
{
  // get the neutronics PDF user object
  FEProblem & from_problem = _multi_app->problem();
  const NeutronicsPDF & neutronics_pdf = from_problem.getUserObject<NeutronicsPDF>(_neutronics_pdf_name);

  // loop over all sub apps and copy over the neutronics data
  for (unsigned int i = 0; i < _multi_app->numGlobalApps(); ++i)
    if (_multi_app->hasLocalApp(i))
    {
      MultiIndex<Real> pdf_data = neutronics_pdf.getPDF(i);

      for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
      {
        PKAGeneratorNeutronicsBase & pka_uo = const_cast<PKAGeneratorNeutronicsBase &>(_multi_app->appProblem(i).getUserObject<PKAGeneratorNeutronicsBase>(_pka_generator_name, tid));
        pka_uo.setPDF(pdf_data);
      }
    }
}
