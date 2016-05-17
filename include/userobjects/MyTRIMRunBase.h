/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MYTRIMRUNBASE_H
#define MYTRIMRUNBASE_H

#include "GeneralUserObject.h"
#include "MyTRIMRasterizer.h"
#include "ThreadedRecoilElementAveragedLoop.h"

#include <map>
#include <vector>

class MyTRIMRunBase;
class MooseMesh;

template<>
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

protected:
  /// Rasterizer object to provide the material data
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements in the TRIM simulation
  const unsigned int _nvars;

  /// number of primary knock-on atoms (PKA) to simulate
  const std::vector<MyTRIM_NS::IonBase> & _pka_list;

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;
};

#endif //MYTRIMRUNBASE_H
