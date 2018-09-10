/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "ThreadedRecoilDiracSourceLoop.h"

ThreadedRecoilDiracSourceLoop::ThreadedRecoilDiracSourceLoop(const MyTRIMRasterizer & rasterizer,
                                                             const MooseMesh & mesh)
  : ThreadedRecoilLoopBase(rasterizer, mesh)
{
}

// Splitting Constructor
ThreadedRecoilDiracSourceLoop::ThreadedRecoilDiracSourceLoop(
    const ThreadedRecoilDiracSourceLoop & x, Threads::split /*split*/)
  : ThreadedRecoilLoopBase(x, Threads::split())
{
}

void
ThreadedRecoilDiracSourceLoop::join(const ThreadedRecoilDiracSourceLoop & rl)
{
  _result_list.insert(_result_list.end(), rl._result_list.begin(), rl._result_list.end());
}

void
ThreadedRecoilDiracSourceLoop::addDefectToResult(const Point & p,
                                                 unsigned int var,
                                                 Real weight,
                                                 ThreadedRecoilDiracSourceLoop::DefectType type)
{
  // TODO: if weight != 1, we need to add a fractional # of results
  if (weight != 1)
    mooseError("Weight != 1 is currently not supported in ThreadedRecoilDiracSourceLoop. This "
               "usually occurs when setting analytical_cutoff != 0.");

  const Elem * elem = (*_pl)(p);
  if (elem != nullptr && var < _nvars)
    _result_list.push_back(MyTRIMResult(p, var, type, elem->id()));
}

template <>
void
dataStore(std::ostream & stream, ThreadedRecoilDiracSourceLoop::MyTRIMResult & dsl, void * context)
{
  dataStore(stream, dsl._location, context);
  dataStore(stream, dsl._var, context);
  dataStore(stream, dsl._type, context);
  dataStore(stream, dsl._elem_id, context);
}

template <>
void
dataLoad(std::istream & stream, ThreadedRecoilDiracSourceLoop::MyTRIMResult & dsl, void * context)
{
  dataLoad(stream, dsl._location, context);
  dataLoad(stream, dsl._var, context);
  dataLoad(stream, dsl._type, context);
  dataLoad(stream, dsl._elem_id, context);
}
