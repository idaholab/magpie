/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MOOSEMYTRIMCORE_H
#define MOOSEMYTRIMCORE_H

#include "mytrim/trim.h"
#include "MooseTypes.h"

#include <list>

class MooseMyTRIMSample;

/**
 * MyTRIM simulation core class
 */
class MooseMyTRIMCore : public MyTRIM_NS::TrimBase
{
public:
  MooseMyTRIMCore(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample);

  virtual void vacancyCreation();
  virtual void replacementCollision();

protected:
  // dimension of the mesh
  const unsigned int _dim;
};

#endif //MOOSEMYTRIMCORE_H
