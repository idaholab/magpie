/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMRasterizer.h"
#include "MooseMyTRIMMaterial.h"

#include "mytrim/sample.h"
#include "mytrim/ion.h"
#include "mytrim/element.h"

// Forward declarations
class MooseMesh;

/**
 * MyTRIM sample class that uses PointLocator lookups on a MooseMesh to
 * fetch and dynamically prepare material data from a MyTRIMRasterizer.
 * A new MooseMyTRIMSample class has to be constructed for every timestep the
 * TRIM simulation is run if the mesh has changed.
 */
class MooseMyTRIMSample : public MyTRIM_NS::SampleBase
{
public:
  MooseMyTRIMSample(const MyTRIMRasterizer &, const MooseMesh &, MyTRIM_NS::SimconfType *);

  /// average crossections for current ion
  virtual void averages(const MyTRIM_NS::IonBase * pka);

  /// interface called by MyTRIM to look up material data
  virtual MyTRIM_NS::MaterialBase * lookupMaterial(Point & pos);

  /// get the simulation dimension
  unsigned int getDim() { return _dim; }

protected:
  /// the rasterizer provides average concentrations for each element
  const MyTRIMRasterizer & _rasterizer;

  /// trim simulation parameters
  const MyTRIMRasterizer::TrimParameters & _trim_parameters;

  /// number of elements used in the problem
  unsigned int _nvars;

  /// mesh of the simulation domain
  const MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;

  /// point locator to determine element pointers form locations
  UniquePtr<PointLocatorBase> _pl;

  ///@{ material cache maps
  typedef std::map<const Elem *, MooseMyTRIMMaterial> MaterialsCache;
  MaterialsCache _materials_master_cache;
  std::map<MyTRIM_NS::IonBase, MaterialsCache> _per_pka_materials_cache;
  ///@}

  /// current ion (for on-the fly averaging)
  const MyTRIM_NS::IonBase * _current_ion;

  /// internal state object for the TRIM simulation
  MyTRIM_NS::SimconfType * _simconf;
};

