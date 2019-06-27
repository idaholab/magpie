/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralUserObject.h"
#include "MyTRIMRasterizer.h"

#include <map>
#include <vector>

class MyTRIMRunBase;
class MooseMesh;

template <>
InputParameters validParams<MyTRIMRunBase>();

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library
 */
class MyTRIMRunBase : public GeneralUserObject
{
public:
  MyTRIMRunBase(const InputParameters & parameters);

  // get the number of elements in the TRIM simulation
  unsigned int nVars() const { return _nvars; }

  // get the name of the rasterizer object (this is to add an artificial dependence in an AuxKernel)
  UserObjectName getRasterizerName() const { return getParam<UserObjectName>("rasterizer"); }

  // get const reference to the rasterizer object
  const MyTRIMRasterizer & rasterizer() const { return _rasterizer; }

protected:
  /// Rasterizer object to provide the material data
  const MyTRIMRasterizer & _rasterizer;

  /// trim simulation parameters
  const MyTRIMRasterizer::TrimParameters & _trim_parameters;

  /// number of elements in the TRIM simulation
  const unsigned int _nvars;

  /// number of primary knock-on atoms (PKA) to simulate
  const std::vector<MyTRIM_NS::IonBase> & _pka_list;

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;
};

