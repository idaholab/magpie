/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MagpieTestApp.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("Magpie");

// Begin the main program.
int
main(int argc, char * argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // print ASCII art application logo
  MagpieTestApp::printLogo();

  // Register this application's MooseApp and any it depends on
  MagpieTestApp::registerApps();

  // This creates dynamic memory that we're responsible for deleting
  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("MagpieTestApp", argc, argv);

  std::shared_ptr<MagpieTestApp> btapp = std::dynamic_pointer_cast<MagpieTestApp>(app);
  if (!btapp)
    mooseError("Did not create Magpie");

  app->setErrorOverridden();

  // Execute the application
  app->run();

  return 0;
}
