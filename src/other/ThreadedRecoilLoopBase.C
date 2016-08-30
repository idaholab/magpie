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
  {
    // seed the RNG with the seed assigned to the primary knock on atom (for parallel reproducibility)
    _simconf.seed(pka._seed);

    // push primary knock on atom onto the recoil queue
    recoils.push(new MyTRIM_NS::IonBase(pka));

    MyTRIM_NS::IonBase * recoil;
    while (!recoils.empty())
    {
      recoil = recoils.front();
      recoils.pop();
      sample.averages(recoil);

      // project into xy plane
      if (_dim == 2)
      {
        recoil->_pos(2) = 0.0;
        recoil->_dir(2) = 0.0;
      }

      // follow this ion's trajectory and store recoils
      // Moose::out << "PKA " << ::round(recoil->_E) << "eV (" << recoil->_Ef << "eV) at " << recoil->_pos(0) << ' ' << recoil->_pos(1) << ' ' << recoil->_pos(2) << ' ' << vac.size() << '\n';
      TRIM.trim(recoil, recoils);

      // store interstitials
      if (recoil->tag >= 0 && recoil->state == MyTRIM_NS::IonBase::INTERSTITIAL)
      {
        // locate element the interstitial is deposited in
        addInterstitialToResult(_rasterizer.periodicPoint(recoil->_pos), recoil->tag);
      }

      // store vacancies
      for (unsigned int i = 0; i < vac.size(); ++i)
        addVacancyToResult(_rasterizer.periodicPoint(vac[i].first), vac[i].second);
      vac.clear();

      // done with this recoil
      delete recoil;
    }
  }
}
