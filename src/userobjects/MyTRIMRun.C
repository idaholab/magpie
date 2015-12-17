/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMRun.h"
#include "MooseMesh.h"
#include "MooseMyTRIMSample.h"
#include "mytrim/trim.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

#include <queue>

template<>
InputParameters validParams<MyTRIMRun>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<UserObjectName>("rasterizer", "MyTRIMRasterizer object to provide material data");
  MultiMooseEnum setup_options(SetupInterface::getExecuteOptions());
  // we run this object once a timestep
  setup_options = "timestep_begin";
  params.set<MultiMooseEnum>("execute_on") = setup_options;
  return params;
}

MyTRIMRun::MyTRIMRun(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _nvars(_rasterizer.nVars()),
    _mesh(_subproblem.mesh()),
    _zero(_nvars, std::pair<Real, Real>(0.0, 0.0))
{
  if (_mesh.isParallelMesh())
    mooseError("MyTRIM runs currently require serial meshes.");

  if (_app.n_processors() > 1)
    mooseError("Parallel communication is not yet implemented. Waiting on libmesh/#748.");
}

void
MyTRIMRun::execute()
{
  // bail out early if no run is requested for this timestep
  if (!_rasterizer.executeThisTimestep())
    return;

  // refresh point locator in case the mesh has changed
  _pl = _mesh.getMesh().sub_point_locator();

  // create a new sample class to bridge the MOOSE mesh and the MyTRIM domain
  MooseMyTRIMSample sample(_rasterizer, _mesh, &_simconf);

  // create a FIFO for recoils
  std::queue<MyTRIM_NS::ionBase *> recoils;

  // use the vacancy mapping TRIM module
  MyTRIM_NS::trimBase TRIM(&_simconf, &sample);

  // create an ion
  MyTRIM_NS::ionBase * pka = new MyTRIM_NS::ionBase;
  pka->gen = 0;  // generation (0 = PKA)
  pka->tag = -1; // tag holds the element type
  pka->z1 = 20;
  pka->m1 = 40;
  pka->e  = 1000;

  pka->dir[0] = 0.0;
  pka->dir[1] = 1.0;
  pka->dir[2] = 0.0;

  pka->pos[0] = 0;
  pka->pos[1] = 0.01;
  pka->pos[2] = 0;

  pka->set_ef();
  recoils.push( pka );

  while (!recoils.empty())
  {
    pka = recoils.front();
    recoils.pop();
    sample.averages(pka);

    // follow this ion's trajectory and store recoils
    TRIM.trim(pka, recoils);

    // store results
    if (pka->tag >= 0)
    {
      // locate element the interstitial is deposited in
      Point p(pka->pos[0], pka->pos[1], pka->pos[2]);
      const Elem * elem = (*_pl)(p);
      if (elem != NULL)
      {
        // store into _result_map
        MyTRIMResultMap::iterator i = _result_map.find(elem->id());
        if (i == _result_map.end())
          i = _result_map.insert(_result_map.begin(), std::make_pair(elem->id(), MyTRIMResult(_nvars, std::make_pair(0.0, 0.0))));

        // increase the interstitial counter for the tagged element
        i->second[pka->tag].second += 1.0;
      }
    }

    // done with this recoil
    delete pka;
  }
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
