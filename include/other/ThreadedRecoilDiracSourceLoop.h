#ifndef THREADEDRECOILDIRACSOURCELOOP_H
#define THREADEDRECOILDIRACSOURCELOOP_H

#include "ThreadedRecoilLoopBase.h"

/**
 * MyTRIM simulation threaded loop for recoil calculation. Results are stored as
 * point lists for use with DiracKernel objects.
 */
class ThreadedRecoilDiracSourceLoop : public ThreadedRecoilLoopBase
{
public:
  ThreadedRecoilDiracSourceLoop(const MyTRIMRasterizer &, const MooseMesh &);

  /// Splitting constructor
  ThreadedRecoilDiracSourceLoop(const ThreadedRecoilDiracSourceLoop & x, Threads::split split);

  /// thread join method
  virtual void join(const ThreadedRecoilDiracSourceLoop &);

  /**
   * result data map for the TRIM simulation holding defects and their locations
   * for each species in the rasterizer.
   */
  typedef std::pair<std::pair<DefectType, unsigned int>, Point> MyTRIMResult;
  typedef std::vector<MyTRIMResult> MyTRIMResultList;
  const MyTRIMResultList & getResultList() { return _result_list; }

protected:
  /// add an interstitial or vacancy to the result list
  void addDefectToResult(const Point & p, unsigned int var, DefectType type);

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultList _result_list;
};

#endif //THREADEDRECOILDIRACSOURCELOOP_H
