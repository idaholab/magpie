/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMDiracRun.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

#include <queue>

template<>
InputParameters validParams<MyTRIMDiracRun>()
{
  InputParameters params = validParams<MyTRIMRunBase>();
  params.addClassDescription("Run a TRIM binary collision Monte Carlo simulation across the entire sample and gather the results for use with a source DiracKernel.");
  return params;
}

MyTRIMDiracRun::MyTRIMDiracRun(const InputParameters & parameters) :
    MyTRIMRunBase(parameters)
{
}

void
MyTRIMDiracRun::execute()
{
  // bail out early if no run is requested for this timestep
  if (!_rasterizer.executeThisTimestep())
    return;

  // trigger the creation of a master PointLocatior outside the threaded region
  _mesh.getPointLocator();

  // build thread loop functor
  ThreadedRecoilDiracSourceLoop rl(_rasterizer, _mesh);

  // output the number of recoils being launched
  _console << "\nMyTRIM: Running " << _pka_list.size() << " recoils." << std::endl;

  // run threads
  Moose::perf_log.push("MyTRIMRecoilLoop", "Solve");
  Threads::parallel_reduce(PKARange(_pka_list.begin(), _pka_list.end()), rl);
  Moose::perf_log.pop("MyTRIMRecoilLoop", "Solve");

  // fetch the joined data from thread 0
  _result_list = rl.getResultList();
}

void
MyTRIMDiracRun::finalize()
{
  // create a one send buffer for use with the libMesh packed range routines
  std::vector<std::string> send_buffers(1);

  // create byte buffers for the streams received from all processors
  std::vector<std::string> recv_buffers;
  recv_buffers.reserve(_app.n_processors());

  // pack the comples datastructures into the string stream
  serialize(send_buffers[0]);

  // broadcast serialized data to and receive from all processors
  _communicator.allgather_packed_range((void *)(NULL), send_buffers.begin(), send_buffers.end(),
                                       std::back_inserter(recv_buffers));

  // unpack the received data and merge it into the local data structures
  deserialize(recv_buffers);
}

const MyTRIMDiracRun::MyTRIMResultList &
MyTRIMDiracRun::result() const
{
  return _result_list;
}

void
MyTRIMDiracRun::serialize(std::string & serialized_buffer)
{
  // stream for serializing the _result_map structure to a byte stream
  std::ostringstream oss;
  dataStore(oss, _result_list, this);

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
}

void
MyTRIMDiracRun::deserialize(std::vector<std::string> & serialized_buffers)
{
  mooseAssert(serialized_buffers.size() == _app.n_processors(), "Unexpected size of serialized_buffers: " << serialized_buffers.size());

  // The input string stream used for deserialization
  std::istringstream iss;

  // Loop over all datastructures for all procerssors to perfrom the gather operation
  for (unsigned int rank = 0; rank < serialized_buffers.size(); ++rank)
  {
    // skip the current processor (its data is already in the strutures)
    if (rank == processor_id())
      continue;

    // populate the stream with a new buffer and reset stream state
    iss.str(serialized_buffers[rank]);
    iss.clear();

    // Load the communicated data into temporary structures
    MyTRIMResultList other_result_list;
    dataLoad(iss, other_result_list, this);

    _result_list.insert(_result_list.end(), other_result_list.begin(), other_result_list.end());
  }
}
