#include "ThreadedRecoilLoopBase.h"
#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"

ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    _rasterizer(rasterizer),
    _nvars(_rasterizer.nVars()),
    _mesh(mesh),
    _dim(_mesh.dimension())
{
}

// Splitting Constructor
ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const ThreadedRecoilLoopBase & x, Threads::split /*split*/) :
    _rasterizer(x._rasterizer),
    _nvars(x._nvars),
    _mesh(x._mesh),
    _dim(x._dim)
{
}

void
ThreadedRecoilLoopBase::operator() (const PKARange & pka_list)
{
  // fetch a point locator
  _pl = _mesh.getPointLocator();

  // create a new sample class to bridge the MOOSE mesh and the MyTRIM domain
  MooseMyTRIMSample sample(_rasterizer, _mesh, &_simconf);

  // create a FIFO for recoils
  std::queue<MyTRIM_NS::IonBase *> recoils;

  // create a list for vacancies created
  std::vector<std::pair<Point, unsigned int> > vac;

  // use the vacancy mapping TRIM module
  MooseMyTRIMCore TRIM(&_simconf, &sample, vac);

  // copy the pka list into the recoil queue
  for (auto && pka : pka_list)
    recoils.push(new MyTRIM_NS::IonBase(pka));

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
    for (unsigned int i = 0; i < vac.size(); ++i)
      addVacancyToResult(_rasterizer.periodicPoint(vac[i].first), vac[i].second);
    vac.clear();

    // done with this recoil
    delete pka;
  }
}
