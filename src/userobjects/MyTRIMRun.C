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

  if (_app.n_processors() > 1)
    mooseError("Parallel communication is not yet implemented. Waiting on libmesh/#748.");

  if (_dim < 2 || _dim > 3)
    mooseError("TRIM simulation works in 2D or 3D only.");

  if (_dim == 2 && (_mesh.getMinInDimension(2) < 0.0 || _mesh.getMaxInDimension(2) > 0.0))
    mooseError("Two dimensional meshes must lie in the z=0 plane.");
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
  /*
  // serialize map into a vector of key,value pairs and allgather it
  typedef std::vector<std::pair<dof_id_type, MyTRIMResult> > Serialization;
  Serialization vecdata(_result_map.begin(), _result_map.end());
  _communicator.allgather(vecdata);

  // sum values on the same elements
  _result_map.clear();
  for (Serialization::const_iterator vec_iter = vecdata.begin(); vec_iter != vecdata.end(); ++vec_iter)
  {
    MyTRIMResultMap::iterator i = _result_map.find(vec_iter->first);

    // if no entry in the map was found then set it, otherwise add value
    if (i == _result_map.end())
      _result_map.insert(*vec_iter);
    else
    {
      mooseAssert(i->second.size() == _nvars && vec_iter->second.size() == _nvars, "Inconsistent TRIM result vector sizes across processors.");

      for (unsigned int j = 0; j < _nvars; ++j)
      {
        // sum vacancy and interstitial production
        i->second[j].first += vec_iter->second[j].first;
        i->second[j].second += vec_iter->second[j].second;
      }
    }
  }
  */
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
