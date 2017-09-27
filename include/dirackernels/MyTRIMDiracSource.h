/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMDIRACSOURCE_H
#define MYTRIMDIRACSOURCE_H

#include "DiracKernel.h"

// forward declarations
class MyTRIMDiracRun;
class MyTRIMDiracSource;

template<>
InputParameters validParams<MyTRIMDiracSource>();

class MyTRIMDiracSource : public DiracKernel
{
public:
  MyTRIMDiracSource(const InputParameters & params);

  virtual void addPoints();
  virtual Real computeQpResidual();

protected:
  const MyTRIMDiracRun & _mytrim;
  const unsigned int _ivar;
  const unsigned int _defect;
};

#endif //MYTRIMDIRACSOURCE_H
