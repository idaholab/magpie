/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef THREADEDRECOILLOOPBASE_H
#define THREADEDRECOILLOOPBASE_H

#include "MooseMyTRIMMaterial.h"
#include "MyTRIMRasterizer.h"
#include "ParallelUniqueId.h"
#include "PointListAdaptor.h"

#include "mytrim/simconf.h"

#include "libmesh/point_locator_base.h"
#include "libmesh/nanoflann.hpp"

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

  /// defect type enum (vacancies, interstitials, and incoming and outgoing replacements)
  enum DefectType {  NONE = -1, VACANCY = 0, INTERSTITIAL, REPLACEMENT_IN, REPLACEMENT_OUT, N_DEFECTS };

protected:
  /// add an interstitial or vacancy to the result list
  virtual void addDefectToResult(const Point & p, unsigned int var, DefectType type) = 0;

  /// add deposited energy to the result list
  virtual void addEnergyToResult(const Point & p, Real edep) = 0;

  /// rasterizer to manage the sample data
  const MyTRIMRasterizer & _rasterizer;

  /// trim simulation parameters
  const MyTRIMRasterizer::TrimParameters & _trim_parameters;

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
