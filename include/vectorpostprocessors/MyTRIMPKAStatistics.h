/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMPKASTATISTICS_H
#define MYTRIMPKASTATISTICS_H

#include "GeneralVectorPostprocessor.h"

// forward declarations
class MyTRIMPKAStatistics;
class MyTRIMRasterizer;

template <>
InputParameters validParams<MyTRIMPKAStatistics>();

class MyTRIMPKAStatistics : public GeneralVectorPostprocessor
{
public:
  MyTRIMPKAStatistics(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

protected:
  /// the rasterizer object to pull data from
  const MyTRIMRasterizer & _rasterizer;

  /// which property is to be aggregated;
  enum ValueType
  {
    MASS = 0,
    ZAID
  } _value_type;

  /// map selected PKA property to the number of PKA with this property
  std::map<unsigned int, unsigned int> _count_map;

  /// property value for the bin
  VectorPostprocessorValue & _property;

  /// number of PKA in the bin
  VectorPostprocessorValue & _count;
};

#endif // MYTRIMPKASTATISTICS_H
