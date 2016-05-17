#ifndef THREADEDRECOILELEMENTAVERAGEDLOOP_H
#define THREADEDRECOILELEMENTAVERAGEDLOOP_H

#include "ThreadedRecoilLoopBase.h"

/**
 * MyTRIM simulation threaded loop for recoil calculation. Results are accumulated
 * as element averages.
 */
class ThreadedRecoilElementAveragedLoop : public ThreadedRecoilLoopBase
{
public:
  ThreadedRecoilElementAveragedLoop(const MyTRIMRasterizer &, const MooseMesh &);

  /// Splitting constructor
  ThreadedRecoilElementAveragedLoop(const ThreadedRecoilElementAveragedLoop & x, Threads::split split);

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
  /// add an interstitial or vacancy to the result list
  void addDefectToResult(const Point & p, unsigned int var, DefectType type);

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultMap _result_map;
};

#endif //THREADEDRECOILELEMENTAVERAGEDLOOP_H
