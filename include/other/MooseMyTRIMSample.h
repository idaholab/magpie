#ifndef MOOSEMYTRIMSAMPLE_H
#define MOOSEMYTRIMSAMPLE_H

#include "MyTRIMRasterizer.h"

#include "mytrim/sample.h"
#include "mytrim/material.h"

// Forward declarations
class MooseMesh;

/**
 * MyTRIM sample class that uses PointLocator lookups on a MooseMesh to
 * fetch and dynamically prepare material data from a MyTRIMRasterizer.
 */
class MooseMyTRIMSample : public MyTRIM_NS::sampleBase
{
public:
  MooseMyTRIMSample(const MyTRIMRasterizer &, const MooseMesh &);

  /// run once per MOOSE timestep
  void initialize();

  /// interface called by MyTRIM to look up material data
  virtual MyTRIM_NS::materialBase* lookupMaterial(double * pos);

protected:

  /// the rasterizer provides average concentrations for each element
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements used in the problem
  unsigned int _nvars;

  ///@{ Element data
  const std::vector<Real> & _trim_mass;
  const std::vector<Real> & _trim_charge;
  ///@}

  /// mesh of the simulation domain
  const MooseMesh & _mesh;

  /// point locator to determine element pointers form locations
  UniquePtr<PointLocatorBase> _pl;

  /// material cache map
  typedef std::map<const Elem *, MyTRIM_NS::materialBase> MaterialsCache;
  MaterialsCache _materials_cache;
};

#endif //MOOSEMYTRIMSAMPLE_H
