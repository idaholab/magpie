/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMElementRun.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

#include <queue>

registerMooseObject("MagpieApp", MyTRIMElementRun);

InputParameters
MyTRIMElementRun::validParams()
{
  InputParameters params = MyTRIMRunBase::validParams();
  params.addClassDescription("Run a TRIM binary collision Monte Carlo simulation across the entire "
                             "sample and gather the results for use with an element averaged "
                             "source Kernel.");
  return params;
}

MyTRIMElementRun::MyTRIMElementRun(const InputParameters & parameters)
  : MyTRIMRunBase(parameters),
    _perf_trim(registerTimedSection("trim", 2)),
    _perf_finalize(registerTimedSection("finalize", 2)),
    _zero(_nvars)
{
}

void
MyTRIMElementRun::execute()
{
  _result_map.clear();

  // bail out early if no run is requested for this timestep
  if (!_rasterizer.executeThisTimestep())
    return;

  // trigger the creation of a master PointLocator outside the threaded region
  _mesh.getPointLocator();

  // build thread loop functor
  ThreadedRecoilElementAveragedLoop rl(_rasterizer, _mesh);

  // output the number of recoils being launched
  _console << "\nMyTRIM: Running " << _trim_parameters.scaled_npka << " recoils.\n"
           << "Recoil rate scaling factor is " << _trim_parameters.recoil_rate_scaling << "\n"
           << "Sampled number of recoils: " << _trim_parameters.original_npka << "\n"
           << "Result scaling factor: " << _trim_parameters.result_scaling_factor << std::endl;

  // run threads
  {
    TIME_SECTION(_perf_trim);
    Threads::parallel_reduce(PKARange(_pka_list.begin(), _pka_list.end()), rl);
  }

  // fetch the joined data from thread 0
  _result_map = rl.getResultMap();
}

void
MyTRIMElementRun::finalize()
{
  TIME_SECTION(_perf_finalize);

  if (!_rasterizer.executeThisTimestep())
    return;

  // for single processor runs we do not need to do anything here
  if (_app.n_processors() > 1)
  {
    // create send buffer
    std::string send_buffer;

    // create byte buffers for the streams received from all processors
    std::vector<std::string> recv_buffers;

    // pack the complex datastructures into the string stream
    serialize(send_buffer);

    // broadcast serialized data to and receive from all processors
    _communicator.allgather(send_buffer, recv_buffers);

    // unpack the received data and merge it into the local data structures
    deserialize(recv_buffers);
  }

  // scale result
  for (auto & result : _result_map)
  {
    result.second._energy *= _trim_parameters.result_scaling_factor;

    for (auto k = beginIndex(result.second._defects); k < _nvars; ++k)
      for (std::size_t l = 0; l < ThreadedRecoilDiracSourceLoop::N_DEFECTS; ++l)
        result.second._defects[k][l] *= _trim_parameters.result_scaling_factor;
  }
}

const MyTRIMElementRun::MyTRIMResult &
MyTRIMElementRun::result(const Elem * elem) const
{
  auto i = _result_map.find(elem->id());

  // if no entry in the map was found no collision event happened in the element
  if (i == _result_map.end())
    return _zero;

  return i->second;
}

void
MyTRIMElementRun::serialize(std::string & serialized_buffer)
{
  // stream for serializing the _result_map structure to a byte stream
  std::ostringstream oss;
  dataStore(oss, _result_map, this);

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
}

void
MyTRIMElementRun::deserialize(std::vector<std::string> & serialized_buffers)
{
  mooseAssert(serialized_buffers.size() == _app.n_processors(),
              "Unexpected size of serialized_buffers: " << serialized_buffers.size());

  // The input string stream used for deserialization
  std::istringstream iss;

  // Loop over all datastructures for all procerssors to perfrom the gather operation
  for (unsigned int rank = 0; rank < serialized_buffers.size(); ++rank)
  {
    // skip the current processor (its data is already in the structures)
    if (rank == processor_id())
      continue;

    // populate the stream with a new buffer and reset stream state
    iss.str(serialized_buffers[rank]);
    iss.clear();

    // Load the communicated data into temporary structures
    MyTRIMResultMap other_result_map;
    dataLoad(iss, other_result_map, this);

    // merge the data in with the current processor's data
    for (auto & other_result : other_result_map)
    {
      auto j = _result_map.lower_bound(other_result.first);

      // if no entry in the map was found then set it, otherwise add value
      if (j == _result_map.end() || j->first != other_result.first)
        _result_map.emplace_hint(j, other_result);
      else
      {
        mooseAssert(other_result.second._defects.size() == j->second._defects.size(),
                    "Inconsistent TRIM defect vector sizes across processors.");
        mooseAssert(j->second._defects.size() == _nvars, "Defect vector size must be _nvars.");

        for (unsigned int k = 0; k < _nvars; ++k)
          for (std::size_t l = 0; l < ThreadedRecoilDiracSourceLoop::N_DEFECTS; ++l)
            j->second._defects[k][l] += other_result.second._defects[k][l];

        // sum energy contributions
        j->second._energy += other_result.second._energy;
      }
    }
  }
}
