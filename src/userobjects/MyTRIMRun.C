/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MyTRIMRun.h"
#include "MooseMesh.h"
#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMSample.h"

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
    _periodic(_rasterizer.periodic()),
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

  // refresh point locator in case the mesh has changed
  _pl = _mesh.getMesh().sub_point_locator();

  // create a new sample class to bridge the MOOSE mesh and the MyTRIM domain
  MooseMyTRIMSample sample(_rasterizer, _mesh, &_simconf);

  // create a FIFO for recoils
  std::queue<MyTRIM_NS::IonBase *> recoils;

  // create a list for vacancies created
  std::vector<std::pair<Point, unsigned int> > vac;

  // use the vacancy mapping TRIM module
  MooseMyTRIMCore TRIM(&_simconf, &sample, vac);

  // copy the pka list into the recoil queue
  for (unsigned int i = 0; i < _pka_list.size(); ++i)
    recoils.push(new MyTRIM_NS::IonBase(_pka_list[i]));

  Moose::out << "Running " << _pka_list.size() << " recoils.\n";

  MyTRIM_NS::IonBase * pka;
  while (!recoils.empty())
  {
    pka = recoils.front();
    recoils.pop();
    sample.averages(pka);

    // project into xy plane
    if (_dim == 2)
    {
      pka->_pos(2) = 0.0;
      pka->_dir(2) = 0.0;
    }

    // follow this ion's trajectory and store recoils
    // Moose::out << "PKA " << ::round(pka->_E) << "eV (" << pka->_Ef << "eV) at " << pka->_pos(0) << ' ' << pka->_pos(1) << ' ' << pka->_pos(2) << ' ' << vac.size() << '\n';
    TRIM.trim(pka, recoils);

    // store interstitials
    if (pka->tag >= 0)
    {
      // locate element the interstitial is deposited in
      addInterstitialToResult(_rasterizer.periodicPoint(pka->_pos), pka->tag);
    }

    // store vacancies
    for (unsigned int i = 0; i < vac.size(); ++i)
      addVacancyToResult(_rasterizer.periodicPoint(vac[i].first), vac[i].second);
    vac.clear();

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

void
MyTRIMRun::addDefectToResult(const Point & p, unsigned int var, MyTRIMRun::DefectType type)
{
  const Elem * elem = (*_pl)(p);
  if (elem != NULL && var < _nvars)
  {
    // store into _result_map
    MyTRIMResultMap::iterator i = _result_map.find(elem->id());
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
