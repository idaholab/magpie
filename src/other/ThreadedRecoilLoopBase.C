#include "ThreadedRecoilLoopBase.h"
#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMEnergyDeposition.h"
#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"

ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    _rasterizer(rasterizer),
    _nvars(_rasterizer.nVars()),
    _mesh(mesh),
    _dim(_mesh.dimension()),
    _analytical_cutoff(_rasterizer.analyticalCutoff()),
    _trim_mass(_rasterizer.mass()),
    _trim_charge(_rasterizer.charge())
{
  _simconf.setLengthScale(_rasterizer.lengthScale());
}

// Splitting Constructor
ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const ThreadedRecoilLoopBase & x, Threads::split /*split*/) :
    _rasterizer(x._rasterizer),
    _nvars(x._nvars),
    _mesh(x._mesh),
    _dim(x._dim),
    _analytical_cutoff(_rasterizer.analyticalCutoff()),
    _trim_mass(_rasterizer.mass()),
    _trim_charge(_rasterizer.charge())
{
  _simconf.setLengthScale(_rasterizer.lengthScale());
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
  std::list<std::pair<Point, unsigned int> > vac_list;

  // create a list potentially used for energy deposition
  std::list<std::pair<Point, Real> > edep_list;

  // build the requested TRIM module
  std::unique_ptr<MooseMyTRIMCore> TRIM;
  switch (_rasterizer.trimModule())
  {
    // basic module with interstitial and vacancy generation
    case MyTRIMRasterizer::MYTRIM_CORE:
      TRIM.reset(new MooseMyTRIMCore(&_simconf, &sample, vac_list));
      break;

    // record energy deposited to the lattice
    case MyTRIMRasterizer::MYTRIM_ENERGY_DEPOSITION:
      TRIM.reset(new MooseMyTRIMEnergyDeposition(&_simconf, &sample, vac_list, edep_list));
      break;

    default:
      mooseError("Unknown TRIM module.");
  }

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

      // full recoil or analytical approximation
      if (recoil->_E < _analytical_cutoff)
      {
        const auto pp = _rasterizer.periodicPoint(recoil->_pos);
        const auto elem = (*_pl)(pp);

        mooseDoOnce(mooseWarning("Skipping detailed cascade calculation below cutoff energy."));

        // Parkin, Don M., and C. Alton Coulter. “Total and Net Displacement
        // Functions for Polyatomic Materials.” Journal of Nuclear Materials
        // 101, no. 3 (October 1, 1981): 261–76. doi:10.1016/0022-3115(81)90169-0.

        // get the composition
        const auto & material_data = _rasterizer.material(elem);
        mooseAssert(material_data.size() == _nvars, "Unexpected material data size");

        // add remaining recoil energy
        if (_rasterizer.trimModule() == MyTRIMRasterizer::MYTRIM_ENERGY_DEPOSITION)
          addEnergyToResult(pp, recoil->_E);
      }
      else
      {
        // follow this ion's trajectory and store recoils
        TRIM->trim(recoil, recoils);

        // store interstitials
        if (recoil->_tag >= 0 && recoil->_state == MyTRIM_NS::IonBase::INTERSTITIAL)
        {
          // locate element the interstitial is deposited in
          addInterstitialToResult(_rasterizer.periodicPoint(recoil->_pos), recoil->_tag);
        }

        // store vacancies
        for (auto & vac: vac_list)
          addVacancyToResult(_rasterizer.periodicPoint(vac.first), vac.second);
        vac_list.clear();

        // store energy deposition
        for (auto & edep: edep_list)
          addEnergyToResult(_rasterizer.periodicPoint(edep.first), edep.second);
        edep_list.clear();
      }

      // done with this recoil
      delete recoil;
    }
  }
}
