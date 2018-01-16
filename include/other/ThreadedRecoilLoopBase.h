/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef THREADEDRECOILLOOPBASE_H
#define THREADEDRECOILLOOPBASE_H

#include "ParallelUniqueId.h"
#include "MooseMyTRIMMaterial.h"
#include "PointListAdaptor.h"

#include "mytrim/simconf.h"

#include "libmesh/point_locator_base.h"
#include "libmesh/nanoflann.hpp"

class MyTRIMRasterizer;
class MooseMesh;

typedef StoredRange<std::vector<MyTRIM_NS::IonBase>::const_iterator, MyTRIM_NS::IonBase> PKARange;

/**
 * MyTRIM simulation threaded loop for recoil calculation. Results are stored as
 * point lists for use with DiracKernel objects.
 */
class ThreadedRecoilLoopBase
{
public:
  ThreadedRecoilLoopBase(const MyTRIMRasterizer &, const MooseMesh &);

  /// Splitting constructor
  ThreadedRecoilLoopBase(const ThreadedRecoilLoopBase & x, Threads::split split);

  /// dummy virtual destructor
  virtual ~ThreadedRecoilLoopBase() {}

  /// parens operator with the code that is executed in threads
  void operator() (const PKARange & range);

  using MyTRIMDefectBufferItem = std::pair<Point, unsigned int>;

protected:
  /// defect type enum
  enum DefectType { VACANCY, INTERSTITIAL, NONE };

  /// add an interstitial or vacancy to the result list
  virtual void addDefectToResult(const Point & p, unsigned int var, DefectType type) = 0;

  /// add deposited energy to the result list
  virtual void addEnergyToResult(const Point & p, Real edep) = 0;

  /// rasterizer to manage the sample data
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements in the TRIM simulation
  unsigned int _nvars;

  /// The Mesh we're using
  const MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;

  /// point locator to use
  UniquePtr<PointLocatorBase> _pl;

  /// internal TRIM simulation status object
  MyTRIM_NS::SimconfType _simconf;

  /// ID number of the current thread
  THREAD_ID _tid;

  /// energy cutoff below which recoils are not followed explicitly but effects are calculated analytically
  const Real _analytical_cutoff;

  ///@{ Element data for Kinchin-Pease
  const std::vector<Real> & _trim_mass;
  const std::vector<Real> & _trim_charge;
  ///@}

private:
  ///@{ Buffer vacancies and interstitials from the same cascade for instantaneous recombination
  std::vector<MyTRIMDefectBufferItem> _vacancy_buffer;
  std::vector<MyTRIMDefectBufferItem> _interstitial_buffer;
  ///@}

  /// add an interstitial to the defect buffer
  void addInterstitialToBuffer(const Point & p, unsigned int var);

  /// add a vacancy to the defect buffer
  void addVacancyToBuffer(const Point & p, unsigned int var);

  using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<Real, PointListAdaptor<MyTRIMDefectBufferItem>>,
      PointListAdaptor<MyTRIMDefectBufferItem>,
      LIBMESH_DIM>;
};

#endif //THREADEDRECOILLOOPBASE_H
