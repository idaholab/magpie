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
#include "ThreadedRecoilLoopBase.h"
#include "MyTRIMRasterizer.h"

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
  const MyTRIMRasterizer & _rasterizer;

  /// rasterizer variable index
  const unsigned int _ivar;

  // defect type to select from the result set for insertion as a source term
  ThreadedRecoilLoopBase::DefectType _defect;

  /// Simulation parameters
  const MyTRIMRasterizer::TrimParameters & _trim_parameters;
};

#endif //MYTRIMDIRACSOURCE_H
