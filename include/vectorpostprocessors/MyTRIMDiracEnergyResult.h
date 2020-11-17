/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralVectorPostprocessor.h"

class MyTRIMDiracEnergyResult;
class MyTRIMDiracRun;

template <>
InputParameters validParams<MyTRIMDiracEnergyResult>();

/**
 * Outputs the list of MyTRIM energy deposition events
 * compiled with the the MyTRIMDiracRunner
 */
class MyTRIMDiracEnergyResult : public GeneralVectorPostprocessor
{
public:
  MyTRIMDiracEnergyResult(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

protected:
  const MyTRIMDiracRun & _mytrim;

  ///@{ coordinates of the defects
  VectorPostprocessorValue & _x;
  VectorPostprocessorValue & _y;
  VectorPostprocessorValue & _z;
  VectorPostprocessorValue & _elem_id;
  VectorPostprocessorValue & _energy_deposition;
  ///@}
};
