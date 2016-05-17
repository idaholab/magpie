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
  for (auto && i : _result_map)
  {
    auto j = _result_map.find(i.first);
    if (j == _result_map.end())
      j = _result_map.insert(_result_map.begin(), std::make_pair(i.first, MyTRIMResult(_nvars, std::make_pair(0.0, 0.0))));

    const MyTRIMResult & src = i.second;
    MyTRIMResult & dst = j->second;

    for (unsigned int k = 0; k < _nvars; ++k)
    {
      // accumulate vacancies and interstitials
      dst[k].first += src[k].first;
      dst[k].second += src[k].second;
    }
  }
}

void
ThreadedRecoilElementAveragedLoop::addDefectToResult(const Point & p, unsigned int var, ThreadedRecoilElementAveragedLoop::DefectType type)
{
  const Elem * elem = (*_pl)(p);
  if (elem != nullptr && var < _nvars)
  {
    // store into _result_map
    auto i = _result_map.find(elem->id());
    if (i == _result_map.end())
      i = _result_map.insert(_result_map.begin(), std::make_pair(elem->id(), MyTRIMResult(_nvars, std::make_pair(0.0, 0.0))));

    // increase the interstitial counter for the tagged element
    switch (type)
    {
      case VACANCY:
        i->second[var].first += 1.0;
        break;

      case INTERSTITIAL:
        i->second[var].second += 1.0;
        break;

      default:
        mooseError("Internal error");
    }
  }
}
