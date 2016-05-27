#ifndef THREADEDRECOILDIRACSOURCELOOP_H
#define THREADEDRECOILDIRACSOURCELOOP_H

#include "ThreadedRecoilLoopBase.h"
#include "DataIO.h"

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
  struct MyTRIMResult;
  typedef std::vector<MyTRIMResult> MyTRIMResultList;
  const MyTRIMResultList & getResultList() { return _result_list; }

protected:
  /// add an interstitial or vacancy to the result list
  void addDefectToResult(const Point & p, unsigned int var, DefectType type);

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultList _result_list;
};

struct ThreadedRecoilDiracSourceLoop::MyTRIMResult {
  MyTRIMResult(const Point & location,
               unsigned int var,
               ThreadedRecoilDiracSourceLoop::DefectType type,
               dof_id_type elem_id) :
      _location(location),
      _var(var),
      _type(type),
      _elem_id(elem_id)
  {
  }

  MyTRIMResult() :
      _location(),
      _var(0),
      _type(NONE),
      _elem_id(libMesh::invalid_uint)
  {
  }

  Point _location;
  unsigned int _var;
  DefectType _type;
  dof_id_type _elem_id;
};

template<>
void dataStore(std::ostream &, ThreadedRecoilDiracSourceLoop::MyTRIMResult &, void *);

template<>
void dataLoad(std::istream &, ThreadedRecoilDiracSourceLoop::MyTRIMResult &, void *);

#endif //THREADEDRECOILDIRACSOURCELOOP_H
