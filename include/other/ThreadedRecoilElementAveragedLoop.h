#ifndef MOOSEMYTRIMTHREADEDRECOILLOOP_H
#define MOOSEMYTRIMTHREADEDRECOILLOOP_H

#include "ParallelUniqueId.h"
#include "MooseMyTRIMMaterial.h"
#include "mytrim/simconf.h"

class MyTRIMRasterizer;
class MooseMesh;

typedef StoredRange<std::vector<MyTRIM_NS::IonBase>::const_iterator, MyTRIM_NS::IonBase> PKARange;

/**
 * MyTRIM simulation threaded loop for recoil calculation
 */
class ThreadedRecoilElementAveragedLoop
{
public:
  ThreadedRecoilElementAveragedLoop(const MyTRIMRasterizer &, const MooseMesh &);

  /// Splitting constructor
  ThreadedRecoilElementAveragedLoop(const ThreadedRecoilElementAveragedLoop & x, Threads::split split);

  /// dummy virtual destructor
  virtual ~ThreadedRecoilElementAveragedLoop() {}

  /// parens operator with the code that is executed in threads
  void operator() (const PKARange & range);

  /// thread join method
  virtual void join(const ThreadedRecoilElementAveragedLoop &);

  /**
   * result data map for the TRIM simulation holding interstitial/vacancy pairs
   * for each species in the rasterizer.
   */
  typedef std::vector<std::pair<Real, Real> > MyTRIMResult;
  typedef std::map<dof_id_type, MyTRIMResult> MyTRIMResultMap;

  const MyTRIMResultMap & getResultMap() { return _result_map; }

protected:
  /// defect type enum
  enum DefectType { VACANCY, INTERSTITIAL };

  /// add an interstitial or vacancy to the result list
  void addDefectToResult(const Point & p, unsigned int var, DefectType type);
  /// add an interstitial to the result list
  void addInterstitialToResult(const Point & p, unsigned int var) { addDefectToResult(p, var, INTERSTITIAL); }
  /// add a vacancy to the result list
  void addVacancyToResult(const Point & p, unsigned int var) { addDefectToResult(p, var, VACANCY); }

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultMap _result_map;

  /// rasterizer to manage the sample data
  const MyTRIMRasterizer & _rasterizer;

  /// number of elements in the TRIM simulation
  unsigned int _nvars;

  /// The Mesh we're using
  const MooseMesh & _mesh;

  /// point locator to use
  UniquePtr<PointLocatorBase> _pl;

  /// dimension of the mesh
  const unsigned int _dim;

  /// internal TRIM simulation status object
  MyTRIM_NS::SimconfType _simconf;

  /// ID number of teh current thread
  THREAD_ID _tid;
};

#endif //MOOSEMYTRIMTHREADEDRECOILLOOP_H
