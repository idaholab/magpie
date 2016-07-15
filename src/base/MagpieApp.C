#include "MagpieApp.h"
#include "Moose.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// AuxKernels
#include "MyTRIMDensityAux.h"
#include "MyTRIMElementResultAux.h"
#include "SPPARKSAux.h"

// Kernels
#include "DefectAnnihilation.h"
#include "MyTRIMElementSource.h"

// DiracKernels
#include "MyTRIMDiracSource.h"

// UserObjects
#include "MyTRIMRasterizer.h"
#include "MyTRIMDiracRun.h"
#include "MyTRIMElementRun.h"
#include "PKAConstant.h"
#include "PKAFissionFragmentEmpirical.h"
#if RATTLESNAKE_ENABLED
  #include "RadiationDamageSN.h"
  #include "RadiationDamageFission.h"
#endif
#include "SPPARKSUserObject.h"

template<>
InputParameters validParams<MagpieApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

MagpieApp::MagpieApp(const InputParameters & parameters) :
    MooseApp(parameters)
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

extern "C" void MagpieApp__registerApps() { MagpieApp::registerApps(); }
void
MagpieApp::registerApps()
{
  registerApp(MagpieApp);
}

void
MagpieApp::registerObjects(Factory & factory)
{
  registerAux(MyTRIMDensityAux);
  registerAux(MyTRIMElementResultAux);
  registerAux(SPPARKSAux);

  registerKernel(DefectAnnihilation);
  registerKernel(MyTRIMElementSource);

  registerDiracKernel(MyTRIMDiracSource);

  registerUserObject(MyTRIMRasterizer);
  registerUserObject(MyTRIMDiracRun);
  registerUserObject(MyTRIMElementRun);
  registerUserObject(PKAConstant);
  registerUserObject(PKAFissionFragmentEmpirical);
  #if RATTLESNAKE_ENABLED
    registerUserObject(RadiationDamageSN);
    registerUserObject(RadiationDamageFission);
  #endif
  registerUserObject(SPPARKSUserObject);
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
