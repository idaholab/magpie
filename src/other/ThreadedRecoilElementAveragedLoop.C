#include "ThreadedRecoilElementAveragedLoop.h"
#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"

ThreadedRecoilElementAveragedLoop::ThreadedRecoilElementAveragedLoop(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    _rasterizer(rasterizer),
    _nvars(_rasterizer.nVars()),
    _mesh(mesh),
    _dim(_mesh.dimension())
{
}

// Splitting Constructor
ThreadedRecoilElementAveragedLoop::ThreadedRecoilElementAveragedLoop(const ThreadedRecoilElementAveragedLoop & x, Threads::split /*split*/) :
    _rasterizer(x._rasterizer),
    _nvars(x._nvars),
    _mesh(x._mesh),
    _dim(x._dim)
{
}

void
ThreadedRecoilElementAveragedLoop::operator() (const PKARange & pka_list)
{
  // fetch a point locator
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
  for (auto && i : pka_list)
    recoils.push(new MyTRIM_NS::IonBase(i));

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
    if (pka->tag >= 0 && pka->state == MyTRIM_NS::IonBase::INTERSTITIAL)
    {
      // locate element the interstitial is deposited in
      addInterstitialToResult(_rasterizer.periodicPoint(pka->_pos), pka->tag);
    }

    // store vacancies
    for (auto && v : vac)
      addVacancyToResult(_rasterizer.periodicPoint(v.first), v.second);
    vac.clear();

    // done with this recoil
    delete pka;
  }
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
