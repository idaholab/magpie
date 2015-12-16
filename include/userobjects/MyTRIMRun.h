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
#include "mytrim/material.h"

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

  /**
   * result data map for the TRIM simulation holding interstitial/vacancy pairs
   * for each species in the rasterizer.
   */
  typedef std::vector<std::pair<Real, Real> > MyTRIMResult;
  typedef std::map<dof_id_type, MyTRIMResult> MyTRIMResultMap;

  /// get the TRIM result data
  const MyTRIMResult & result(const Elem *) const;

protected:
  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultMap _result_map;

  /// Rasterizer object to provide the material data
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements in the TRIM simulation
  int _nvars;

  /// variable number to use for minPeriodicDistance calls (i.e. use the periodicity of this variable)
  int _periodic;

  /// point locator to use
  UniquePtr<PointLocatorBase> _pl;

  /// The Mesh we're using
  MooseMesh & _mesh;

private:
  /// zero result to return for elements that have not been touched by the cascades
  MyTRIMResult _zero;
};

#endif //MYTRIMRUN_H
