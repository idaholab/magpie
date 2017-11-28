/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "AppFactory.h"
#include "MagpieApp.h"
#include "Moose.h"
#include "MooseSyntax.h"

// AuxKernels
#include "AtomicDensityAux.h"
#include "MyTRIMDensityAux.h"
#include "MyTRIMElementEnergyAux.h"
#include "MyTRIMElementResultAux.h"
#include "SPPARKSAux.h"

// Kernels
#include "CoupledDefectAnnihilation.h"
#include "DefectAnnihilation.h"
#include "MyTRIMElementHeatSource.h"
#include "MyTRIMElementSource.h"

// DiracKernels
#include "MyTRIMDiracSource.h"

// UserObjects
#include "MyTRIMDiracRun.h"
#include "MyTRIMElementRun.h"
#include "MyTRIMPKAInfo.h"
#include "MyTRIMPKAInConeInfo.h"
#include "MyTRIMRasterizer.h"
#include "NeutronicsSpectrumSamplerFission.h"
#include "NeutronicsSpectrumSamplerSN.h"
#include "PKAConstant.h"
#include "PKAFissionFragmentEmpirical.h"
#include "PKAFissionFragmentNeutronics.h"
#include "PKAFixedPointGenerator.h"
#include "PKAGun.h"
#include "SPPARKSUserObject.h"
#include "ElasticRecoilCrossSectionUserObject.h"
#include "IsotopeRecoilRate.h"

// Transfers
#include "MultiAppNeutronicsSpectrumTransfer.h"

// Meshes
#include "MyTRIMMesh.h"

// VectorPostprocessors
#include "MyTRIMDiracResult.h"
#include "MyTRIMPKAEnergyHistogram.h"
#include "MyTRIMPKAStatistics.h"
#include "PKAList.h"
#include "IsotopeRecoilRateSampler.h"

template <>
InputParameters
validParams<MagpieApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

MagpieApp::MagpieApp(const InputParameters & parameters)
  : MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  MagpieApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  MagpieApp::associateSyntax(_syntax, _action_factory);
}

MagpieApp::~MagpieApp()
{
}

extern "C" void
MagpieApp__registerApps()
{
  MagpieApp::registerApps();
}
void
MagpieApp::registerApps()
{
  registerApp(MagpieApp);
}

void
MagpieApp::registerObjects(Factory & factory)
{
  registerAux(AtomicDensityAux);
  registerAux(MyTRIMDensityAux);
  registerAux(MyTRIMElementEnergyAux);
  registerAux(MyTRIMElementResultAux);
  registerAux(SPPARKSAux);

  registerKernel(CoupledDefectAnnihilation);
  registerKernel(DefectAnnihilation);
  registerKernel(MyTRIMElementHeatSource);
  registerKernel(MyTRIMElementSource);

  registerDiracKernel(MyTRIMDiracSource);

  registerUserObject(MyTRIMRasterizer);
  registerUserObject(MyTRIMDiracRun);
  registerUserObject(MyTRIMElementRun);
  registerUserObject(PKAConstant);
  registerUserObject(PKAFixedPointGenerator);
  registerUserObject(PKAGun);
  registerUserObject(PKAFissionFragmentEmpirical);
  registerUserObject(PKAFissionFragmentNeutronics);
#ifdef RATTLESNAKE_ENABLED
  registerUserObject(NeutronicsSpectrumSamplerSN);
#endif
  registerUserObject(NeutronicsSpectrumSamplerFission);
  registerUserObject(SPPARKSUserObject);
  registerUserObject(MyTRIMPKAInfo);
  registerUserObject(MyTRIMPKAInConeInfo);
  registerUserObject(ElasticRecoilCrossSectionUserObject);
  registerUserObject(IsotopeRecoilRate);

  registerTransfer(MultiAppNeutronicsSpectrumTransfer);

  registerMesh(MyTRIMMesh);

  registerVectorPostprocessor(MyTRIMDiracResult);
  registerVectorPostprocessor(MyTRIMPKAEnergyHistogram);
  registerVectorPostprocessor(MyTRIMPKAStatistics);
  registerVectorPostprocessor(PKAList);
  registerVectorPostprocessor(IsotopeRecoilRateSampler);
}

void
MagpieApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

void
MagpieApp::printLogo()
{
  Moose::out << "\n"
             << "\n             M   M   AA    GGG  PPP   I   EE   "
             << "\n             MM MM  A  A  G     P  P  I  E  E  "
             << "\n             M M M  AAAA  G  G  PPP   I  EE    "
             << "\n             M   M  A  A   GG   P     I   EEE  "
             << "\n"
             << "\n  Mesoscale Atomistic Glue Program for Integrated Execution"
             << "\n"
             << "\n         `;+;`                                                       "
             << "\n       ,#######`                                                     "
             << "\n     ,##########`                                                    "
             << "\n .;;#############                                                    "
             << "\n`####+'+##########                                                   "
             << "\n  `,+#+############:`                                                "
             << "\n      `#################'`                                           "
             << "\n       :###################+                                         "
             << "\n        ######     .`;#######+                                       "
             << "\n         #'##          ;#+++'##`                                     "
             << "\n         ;+#+            ,;;;;;';                                    "
             << "\n          ###'.,:,` `:;;;;;;;;;;;;;                                  "
             << "\n          ##.'+;;;;;;;;;;;;;;;;;;;;;;`                               "
             << "\n          ## ++;;;;;;;;;;+;'';;;;;;;;#+`                             "
             << "\n          ;#`,;;;;;;;;;;;;,.      .''+##+;`                          "
             << "\n           #' `;;;;.::;;'+'+#######+++';;;;;';.                      "
             << "\n           +#              .;+';';;############+###+;.               "
             << "\n           ,#,               `;;'';#######,;+############+;.         "
             << "\n            ;#               ;+##+######'         `.###########+;`   "
             << "\n              '              +#########                    ,.+######`"
             << "\n                             ########`                           .;:,"
             << "\n                            ,######:                                 "
             << "\n                 ;          #######                                  "
             << "\n                  .;       ,######                                   "
             << "\n                    ,;   .`+#####:                                   "
             << "\n                      ';##;;;+###                                    "
             << "\n                      ;'#  `;;#+                                     "
             << "\n                   :;+',     '+                                      "
             << "\n                .;;'        .#                                       "
             << "\n            `;. ;#,        ,##                                       "
             << "\n             `,,:` +`   .::##' ::.                                   "
             << "\n                   `:'   : +#+.                                      "
             << "\n                     '`      `:                                      "
             << "\n";
}
