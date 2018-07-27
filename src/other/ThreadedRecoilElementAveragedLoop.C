/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "ThreadedRecoilElementAveragedLoop.h"

ThreadedRecoilElementAveragedLoop::ThreadedRecoilElementAveragedLoop(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    ThreadedRecoilLoopBase(rasterizer, mesh)
{
}

// Splitting Constructor
ThreadedRecoilElementAveragedLoop::ThreadedRecoilElementAveragedLoop(const ThreadedRecoilElementAveragedLoop & x, Threads::split /*split*/) :
    ThreadedRecoilLoopBase(x, Threads::split())
{
}

void
ThreadedRecoilElementAveragedLoop::join(const ThreadedRecoilElementAveragedLoop & rl)
{
  for (auto && i : rl._result_map)
  {
    // find the result map entry corresponding to teh other thred's position
    // lower_bound gives the first element that os _not_ less than the searched elemen.
    // This is either a match, or the right position for a hinted insertion (fast)
    auto j = _result_map.lower_bound(i.first);
    if (j == _result_map.end() || j->first != i.first)
      j = _result_map.emplace_hint(j, i.first, MyTRIMResult(_nvars));

    const MyTRIMResult & src = i.second;
    MyTRIMResult & dst = j->second;

    mooseAssert(dst._defects.size() == src._defects.size(), "Defect vector sizes inconsistent.");
    mooseAssert(dst._defects.size() == _nvars, "Defect vector size must be _nvars.");

    // accumulate vacancies, interstitials, and replcaements
    for (auto k = beginIndex(dst._defects); k < _nvars; ++k)
      for (std::size_t l = 0; l < N_DEFECTS; ++l)
        dst._defects[k][l] += src._defects[k][l];

    dst._energy += src._energy;
  }
}

void
ThreadedRecoilElementAveragedLoop::addDefectToResult(const Point & p, unsigned int var, ThreadedRecoilElementAveragedLoop::DefectType type)
{
  const Elem * elem = (*_pl)(p);
  if (elem == nullptr || var >= _nvars)
    return;

  // store into _result_map
  auto i = _result_map.lower_bound(elem->id());
  if (i == _result_map.end() || i->first != elem->id())
    i = _result_map.emplace_hint(i, elem->id(), MyTRIMResult(_nvars));

  // check for invalid types
  mooseAssert(type != NONE && type != N_DEFECTS, "Invalid defect type passed to addDefectToResult");

  // increase the interstitial counter for the tagged element
  i->second._defects[var][type] += 1.0;
}

void
ThreadedRecoilElementAveragedLoop::addEnergyToResult(const Point & p, Real edep)
{
  const Elem * elem = (*_pl)(p);
  if (elem == nullptr)
    return;

  // store into _result_map
  auto i = _result_map.lower_bound(elem->id());
  if (i == _result_map.end() || i->first != elem->id())
    i = _result_map.emplace_hint(i, elem->id(), MyTRIMResult(_nvars));

  i->second._energy += edep;
}

template<>
void
dataStore(std::ostream & stream, ThreadedRecoilElementAveragedLoop::MyTRIMResult & eal, void * context)
{
  dataStore(stream, eal._defects, context);
  dataStore(stream, eal._energy, context);
}

template<>
void
dataLoad(std::istream & stream, ThreadedRecoilElementAveragedLoop::MyTRIMResult & eal, void * context)
{
  dataLoad(stream, eal._defects, context);
  dataLoad(stream, eal._energy, context);
}
