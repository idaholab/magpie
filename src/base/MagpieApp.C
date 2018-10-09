/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MagpieApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<MagpieApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

// When using the new Registry system, this line is required so that
// dependent apps know about the MagpieApp label.
registerKnownLabel("MagpieApp");

MagpieApp::MagpieApp(const InputParameters & parameters) : MooseApp(parameters)
{
  srand(processor_id());
  MagpieApp::registerAll(_factory, _action_factory, _syntax);
}

MagpieApp::~MagpieApp() {}

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

// External entry point for object registration
extern "C" void
MagpieApp__registerAll(Factory & factory, ActionFactory & action_factory, Syntax & syntax)
{
  MagpieApp::registerAll(factory, action_factory, syntax);
}
void
MagpieApp::registerAll(Factory & factory, ActionFactory & action_factory, Syntax & syntax)
{
  Registry::registerObjectsTo(factory, {"MagpieApp"});
  Registry::registerActionsTo(action_factory, {"MagpieApp"});
  MagpieApp::associateSyntax(syntax, action_factory);

  Moose::registerAll(factory, action_factory, syntax);
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
