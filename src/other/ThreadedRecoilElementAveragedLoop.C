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
    auto j = _result_map.find(i.first);
    if (j == _result_map.end())
      j = _result_map.insert(_result_map.begin(), std::make_pair(i.first, MyTRIMResult(_nvars)));

    const MyTRIMResult & src = i.second;
    MyTRIMResult & dst = j->second;

    mooseAssert(dst._defects.size() == src._defects.size(), "Defect vector sizes inconsistent.");
    mooseAssert(dst._defects.size() == _nvars, "Defect vector size must be _nvars.");

    for (unsigned int k = 0; k < _nvars; ++k)
    {
      // accumulate vacancies and interstitials
      dst._defects[k]._vacancies += src._defects[k]._vacancies;
      dst._defects[k]._interstitials += src._defects[k]._interstitials;
    }

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
  auto i = _result_map.find(elem->id());
  if (i == _result_map.end())
    i = _result_map.insert(_result_map.begin(), std::make_pair(elem->id(), MyTRIMResult(_nvars)));

  // increase the interstitial counter for the tagged element
  switch (type)
  {
    case VACANCY:
      i->second._defects[var]._vacancies += 1.0;
      break;

    case INTERSTITIAL:
      i->second._defects[var]._interstitials += 1.0;
      break;

    default:
      mooseError("Internal error");
  }
}

void
ThreadedRecoilElementAveragedLoop::addEnergyToResult(const Point & p, Real edep)
{
  const Elem * elem = (*_pl)(p);
  if (elem == nullptr)
    return;

  // store into _result_map
  auto i = _result_map.find(elem->id());
  if (i == _result_map.end())
    i = _result_map.insert(_result_map.begin(), std::make_pair(elem->id(), MyTRIMResult(_nvars)));

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
