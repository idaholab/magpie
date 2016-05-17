#ifndef THREADEDRECOILLOOPBASE_H
#define THREADEDRECOILLOOPBASE_H

#include "ParallelUniqueId.h"
#include "MooseMyTRIMMaterial.h"
#include "mytrim/simconf.h"
#include "libmesh/point_locator_base.h"

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

protected:
  /// defect type enum
  enum DefectType { VACANCY, INTERSTITIAL };

  /// add an interstitial or vacancy to the result list
  virtual void addDefectToResult(const Point & p, unsigned int var, DefectType type) = 0;

  /// add an interstitial to the result list
  void addInterstitialToResult(const Point & p, unsigned int var) { addDefectToResult(p, var, INTERSTITIAL); }
  /// add a vacancy to the result list
  void addVacancyToResult(const Point & p, unsigned int var) { addDefectToResult(p, var, VACANCY); }

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
};

#endif //THREADEDRECOILLOOPBASE_H
