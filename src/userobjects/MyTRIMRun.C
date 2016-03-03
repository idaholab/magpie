/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMRun.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

#include <queue>

template<>
InputParameters validParams<MyTRIMRun>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription("Run a TRIM binary collision Monte Carlo simulation across the entire sample");
  params.addRequiredParam<UserObjectName>("rasterizer", "MyTRIMRasterizer object to provide material data");

  // we run this object once a timestep
  params.set<MultiMooseEnum>("execute_on") = "timestep_begin";
  params.suppressParameter<MultiMooseEnum>("execute_on");

  return params;
}

MyTRIMRun::MyTRIMRun(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _nvars(_rasterizer.nVars()),
    _pka_list(_rasterizer.getPKAList()),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension()),
    _zero(_nvars, std::pair<Real, Real>(0.0, 0.0))
{
  if (_mesh.isParallelMesh())
    mooseError("MyTRIM runs currently require serial meshes.");

  if (_dim < 2 || _dim > 3)
    mooseError("TRIM simulation works in 2D or 3D only.");

  if (_dim == 2)
  {
    // make sure all nodes lie in the xy-plane
    const MeshBase::node_iterator nd_end = _mesh.getMesh().nodes_end();
    for (MeshBase::node_iterator nd = _mesh.getMesh().nodes_begin(); nd != nd_end; ++nd)
      if ((**nd)(2) != 0.0)
        mooseError("Two dimensional meshes must lie in the z=0 plane.");
  }
}

void
MyTRIMRun::execute()
{
  // bail out early if no run is requested for this timestep
  if (!_rasterizer.executeThisTimestep())
    return;

  // trigger the creation of a master PointLocatior outside the threaded region
  _mesh.getMesh().sub_point_locator();

  // build thread loop functor
  MooseMyTRIMThreadedRecoilLoop rl(_rasterizer, _mesh);

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
MyTRIMRun::finalize()
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

const MyTRIMRun::MyTRIMResult &
MyTRIMRun::result(const Elem * elem) const
{
  MyTRIMResultMap::const_iterator i = _result_map.find(elem->id());

  // if no entry in the map was found no collision event happened in the element
  if (i == _result_map.end())
    return _zero;

  return i->second;
}

void
MyTRIMRun::serialize(std::string & serialized_buffer)
{
  // stream for serializing the _result_map structure to a byte stream
  std::ostringstream oss;
  dataStore(oss, _result_map, this);

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
}

void
MyTRIMRun::deserialize(std::vector<std::string> & serialized_buffers)
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
    MyTRIMResultMap other_result_map;
    dataLoad(iss, other_result_map, this);

    // merge the data in with the current processor's data
    for (MyTRIMResultMap::const_iterator i = other_result_map.begin(); i != other_result_map.end(); ++i)
    {
      MyTRIMResultMap::iterator j = _result_map.find(i->first);

      // if no entry in the map was found then set it, otherwise add value
      if (j == _result_map.end())
        _result_map.insert(*i);
      else
      {
        mooseAssert(i->second.size() == _nvars && j->second.size() == _nvars, "Inconsistent TRIM result vector sizes across processors.");

        for (unsigned int k = 0; k < _nvars; ++k)
        {
          // sum vacancy and interstitial production
          j->second[k].first += i->second[k].first;
          j->second[k].second += i->second[k].second;
        }
      }
    }
  }
}
