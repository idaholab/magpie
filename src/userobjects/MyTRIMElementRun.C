#include "MyTRIMElementRun.h"
#include "MagpieParallel.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

#include <queue>

template<>
InputParameters validParams<MyTRIMElementRun>()
{
  InputParameters params = validParams<MyTRIMRunBase>();
  params.addClassDescription("Run a TRIM binary collision Monte Carlo simulation across the entire sample and gather the results for use with an element averaged source Kernel.");
  return params;
}

MyTRIMElementRun::MyTRIMElementRun(const InputParameters & parameters) :
    MyTRIMRunBase(parameters),
    _zero(_nvars, std::pair<Real, Real>(0.0, 0.0))
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
  _console << "\nMyTRIM: Running " << _pka_list.size() << " recoils." << std::endl;

  // run threads
  Moose::perf_log.push("MyTRIMRecoilLoop", "Solve");
  Threads::parallel_reduce(PKARange(_pka_list.begin(), _pka_list.end()), rl);

  Moose::perf_log.pop("MyTRIMRecoilLoop", "Solve");

  // fetch the joined data from thread 0
  _result_map = rl.getResultMap();
}

void
MyTRIMElementRun::finalize()
{
  // for single processor runs we do not need to do anything here
  if (!_rasterizer.executeThisTimestep() || _app.n_processors() == 1)
    return;

  // create send buffer
  std::string send_buffer;

  // create byte buffers for the streams received from all processors
  std::vector<std::string> recv_buffers;

  // pack the complex datastructures into the string stream
  serialize(send_buffer);

  // broadcast serialized data to and receive from all processors
  MagpieUtils::allgatherStringBuffers(_communicator, send_buffer, recv_buffers);

  // unpack the received data and merge it into the local data structures
  deserialize(recv_buffers);
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
  mooseAssert(serialized_buffers.size() == _app.n_processors(), "Unexpected size of serialized_buffers: " << serialized_buffers.size());

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
      auto j = _result_map.find(other_result.first);

      // if no entry in the map was found then set it, otherwise add value
      if (j == _result_map.end())
        _result_map.insert(other_result);
      else
      {
        mooseAssert(other_result.second.size() == _nvars && j->second.size() == _nvars, "Inconsistent TRIM result vector sizes across processors.");

        for (unsigned int k = 0; k < _nvars; ++k)
        {
          // sum vacancy and interstitial production
          j->second[k].first += other_result.second[k].first;
          j->second[k].second += other_result.second[k].second;
        }
      }
    }
  }
}
