/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MagpieTestApp.h"

#include "AppFactory.h"

InputParameters
MagpieTestApp::validParams()
{
  InputParameters params = MagpieApp::validParams();
  return params;
}

MagpieTestApp::MagpieTestApp(InputParameters parameters) : MagpieApp(parameters)
{
  MagpieTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

MagpieTestApp::~MagpieTestApp() {}

// External entry point for dynamic application loading
extern "C" void
MagpieTestApp__registerApps()
{
  MagpieTestApp::registerApps();
}

void
MagpieTestApp::registerApps()
{
  MagpieApp::registerApps();
  registerApp(MagpieTestApp);
}

// External entry point for object registration
extern "C" void
MagpieApp__registerAll(Factory & factory, ActionFactory & action_factory, Syntax & syntax)
{
  MagpieTestApp::registerAll(factory, action_factory, syntax);
}

void
MagpieTestApp::registerAll(Factory & factory,
                           ActionFactory & action_factory,
                           Syntax & syntax,
                           bool use_test_objs)
{
  MagpieApp::registerAll(factory, action_factory, syntax);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(factory, {"MagpieTestApp"});
    Registry::registerActionsTo(action_factory, {"MagpieTestApp"});
    // register MagpieTestApp syntax here if any gets added
  }
}
