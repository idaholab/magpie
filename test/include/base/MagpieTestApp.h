/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MagpieApp.h"

class MagpieTestApp : public MagpieApp
{
public:
  MagpieTestApp(InputParameters parameters);
  static InputParameters validParams();

  virtual ~MagpieTestApp();
  static void registerApps();
  static void registerAll(Factory & factory,
                          ActionFactory & action_factory,
                          Syntax & syntax,
                          bool use_test_objs = false);
};
