/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralVectorPostprocessor.h"

// forward declarations
class MyTRIMPKAEnergyHistogram;
class MyTRIMRasterizer;

template <>
InputParameters validParams<MyTRIMPKAEnergyHistogram>();

class MyTRIMPKAEnergyHistogram : public GeneralVectorPostprocessor
{
public:
  MyTRIMPKAEnergyHistogram(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

protected:
  /// the rasterizer object to pull data from
  const MyTRIMRasterizer & _rasterizer;

  /// number of energy channels
  const unsigned int _nchannels;

  /// bin width
  Real _deltaE;

  /// energy mid point of the channel
  VectorPostprocessorValue & _channel_center;

  /// number of PKA in the channel
  VectorPostprocessorValue & _count;
};

