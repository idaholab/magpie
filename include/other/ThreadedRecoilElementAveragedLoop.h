/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ThreadedRecoilLoopBase.h"
#include "DataIO.h"

/**
 * MyTRIM simulation threaded loop for recoil calculation. Results are accumulated
 * as element averages.
 */
class ThreadedRecoilElementAveragedLoop : public ThreadedRecoilLoopBase
{
public:
  ThreadedRecoilElementAveragedLoop(const MyTRIMRasterizer &, const MooseMesh &);

  /// Splitting constructor
  ThreadedRecoilElementAveragedLoop(const ThreadedRecoilElementAveragedLoop & x,
                                    Threads::split split);

  /// thread join method
  virtual void join(const ThreadedRecoilElementAveragedLoop &);

  /**
   * result data map for the TRIM simulation holding interstitial/vacancy pairs
   * for each species in the rasterizer.
   */
  struct MyTRIMResult;
  typedef std::map<dof_id_type, MyTRIMResult> MyTRIMResultMap;

  const MyTRIMResultMap & getResultMap() { return _result_map; }

protected:
  /// add an interstitial or vacancy to the result list
  void addDefectToResult(const Point & p, unsigned int var, Real weight, DefectType type);

  /// add deposited energy to the result list
  void addEnergyToResult(const Point & p, Real edep);

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultMap _result_map;
};

struct ThreadedRecoilElementAveragedLoop::MyTRIMResult
{
  MyTRIMResult(unsigned int nvars) : _defects(nvars), _energy(0.0) {}
  MyTRIMResult() : _defects(), _energy(0.0) {}

  // numbers of point defects per chemical element
  using Defect = std::array<Real, 4>;
  std::vector<Defect> _defects;

  /// this will hold the matrix of replacement collisions
  // std::vector<Real> _replacements;

  /// deposited energy per element
  Real _energy;
};

template <>
void dataStore(std::ostream &, ThreadedRecoilElementAveragedLoop::MyTRIMResult &, void *);

template <>
void dataLoad(std::istream &, ThreadedRecoilElementAveragedLoop::MyTRIMResult &, void *);
