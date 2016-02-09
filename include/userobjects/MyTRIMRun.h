/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MYTRIMRUN_H
#define MYTRIMRUN_H

#include "GeneralUserObject.h"
#include "MyTRIMRasterizer.h"
#include "MooseMyTRIMThreadedRecoilLoop.h"

#include <map>
#include <vector>

class MyTRIMRun;
class MooseMesh;

template<>
InputParameters validParams<MyTRIMRun>();

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library
 */
class MyTRIMRun : public GeneralUserObject
{
public:
  MyTRIMRun(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute();
  virtual void finalize();

  /// @{ shorthand typedefs
  typedef MooseMyTRIMThreadedRecoilLoop::MyTRIMResult MyTRIMResult;
  typedef MooseMyTRIMThreadedRecoilLoop::MyTRIMResultMap MyTRIMResultMap;
  /// @}

  /// get the TRIM result data
  const MyTRIMResult & result(const Elem *) const;

  // get the number of elements in the TRIM simulation
  unsigned int nVars() const { return _nvars; }

protected:
  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultMap _result_map;

  /// Rasterizer object to provide the material data
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements in the TRIM simulation
  int _nvars;

  /// number of primary knock-on atoms (PKA) to simulate
  const std::vector<MyTRIM_NS::IonBase> & _pka_list;

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;

private:
  /// zero result to return for elements that have not been touched by the cascades
  MyTRIMResult _zero;
};

#endif //MYTRIMRUN_H
