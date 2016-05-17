#include "ThreadedRecoilDiracSourceLoop.h"

ThreadedRecoilDiracSourceLoop::ThreadedRecoilDiracSourceLoop(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    ThreadedRecoilLoopBase(rasterizer, mesh)
{
}

// Splitting Constructor
ThreadedRecoilDiracSourceLoop::ThreadedRecoilDiracSourceLoop(const ThreadedRecoilDiracSourceLoop & x, Threads::split /*split*/) :
    ThreadedRecoilLoopBase(x, Threads::split())
{
}

void
ThreadedRecoilDiracSourceLoop::join(const ThreadedRecoilDiracSourceLoop & rl)
{
  _result_list.insert(_result_list.end(), rl._result_list.begin(), rl._result_list.end());
}

void
ThreadedRecoilDiracSourceLoop::addDefectToResult(const Point & p, unsigned int var, ThreadedRecoilDiracSourceLoop::DefectType type)
{
  _result_list.push_back(std::make_pair(std::make_pair(type, var), p));
}
