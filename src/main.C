#include "MagpieApp.h"
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
  MagpieApp::printLogo();

  // Register this application's MooseApp and any it depends on
  MagpieApp::registerApps();

  // This creates dynamic memory that we're responsible for deleting
  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("MagpieApp", argc, argv);

  // Execute the application
  app->run();

  return 0;
}
