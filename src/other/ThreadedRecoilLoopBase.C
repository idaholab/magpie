/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "ThreadedRecoilLoopBase.h"
#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMEnergyDeposition.h"
#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"

// Specialization for PointListAdaptor<MyTRIMDefectBufferItem>
template <>
inline const Point &
PointListAdaptor<ThreadedRecoilLoopBase::MyTRIMDefectBufferItem>::getPoint(const ThreadedRecoilLoopBase::MyTRIMDefectBufferItem & item) const { return item.first; }

ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    _rasterizer(rasterizer),
    _trim_parameters(_rasterizer.getTrimParameters()),
    _nvars(_trim_parameters.nVars()),
    _mesh(mesh),
    _dim(_mesh.dimension())
{
  _simconf.setLengthScale(_trim_parameters.length_scale);
}

// Splitting Constructor
ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const ThreadedRecoilLoopBase & x, Threads::split /*split*/) :
    _rasterizer(x._rasterizer),
    _trim_parameters(x._trim_parameters),
    _nvars(x._nvars),
    _mesh(x._mesh),
    _dim(x._dim)
{
  _simconf.setLengthScale(_trim_parameters.length_scale);
}

void
ThreadedRecoilLoopBase::operator() (const PKARange & pka_list)
{
  // fetch a point locator
  _pl = _mesh.getPointLocator();

  // permit querying points that are potentially outside the mesh
  _pl->enable_out_of_mesh_mode();

  // create a new sample class to bridge the MOOSE mesh and the MyTRIM domain
  MooseMyTRIMSample sample(_rasterizer, _mesh, &_simconf);

  // create a FIFO for recoils
  std::queue<MyTRIM_NS::IonBase *> recoils;

  // create a list for vacancies created
  std::list<MyTRIMDefectBufferItem> vac_list;

  // create a list potentially used for energy deposition
  std::list<std::pair<Point, Real> > edep_list;

  // build the requested TRIM module
  std::unique_ptr<MooseMyTRIMCore> TRIM;
  switch (_trim_parameters.trim_module)
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

  // for instantaneous recombination
  const unsigned int requested_results = 10;
  std::vector<std::size_t> return_index(requested_results);
  std::vector<Real> return_dist_sqr(requested_results);

  // copy the pka list into the recoil queue
  for (auto && pka : pka_list)
  {
    // seed the RNG with the seed assigned to the primary knock on atom (for parallel reproducibility)
    _simconf.seed(pka._seed);

    // push primary knock on atom onto the recoil queue
    recoils.push(new MyTRIM_NS::IonBase(pka));

    // clear cascade recombination buffers
    _interstitial_buffer.clear();
    _vacancy_buffer.clear();

    MyTRIM_NS::IonBase * recoil;
    while (!recoils.empty())
    {
      recoil = recoils.front();
      recoils.pop();
      sample.averages(recoil);

      // project into xy plane
      if (_dim == 2)
        recoil->_pos(2) = 0.0;

      // full recoil or analytical approximation
      if (recoil->_E < _trim_parameters.analytical_cutoff)
      {
        const auto pp = _rasterizer.periodicPoint(recoil->_pos);
        // const auto elem = (*_pl)(pp);

        mooseDoOnce(mooseWarning("Skipping detailed cascade calculation below cutoff energy."));

        // Parkin, Don M., and C. Alton Coulter. “Total and Net Displacement
        // Functions for Polyatomic Materials.” Journal of Nuclear Materials
        // 101, no. 3 (October 1, 1981): 261–76. doi:10.1016/0022-3115(81)90169-0.

        // get the composition
        // const auto & material_data = _rasterizer.material(elem);
        // mooseAssert(material_data.size() == _nvars, "Unexpected material data size");

        // add remaining recoil energy
        if (_trim_parameters.trim_module == MyTRIMRasterizer::MYTRIM_ENERGY_DEPOSITION)
          addEnergyToResult(pp, recoil->_E);
      }
      else
      {
        // follow this ion's trajectory and store recoils
        TRIM->trim(recoil, recoils);

        // store interstitials
        if (recoil->_tag >= 0 && recoil->_state == MyTRIM_NS::IonBase::INTERSTITIAL)
          _interstitial_buffer.push_back(std::make_pair(recoil->_pos, recoil->_tag));

        // store vacancies
        for (auto & vac: vac_list)
          _vacancy_buffer.push_back(vac);
        vac_list.clear();

        // store energy deposition
        for (auto & edep: edep_list)
          addEnergyToResult(_rasterizer.periodicPoint(edep.first), edep.second);
        edep_list.clear();
      }

      // done with this recoil
      delete recoil;
    }

    //
    // Process instantaneous recombination of this PKA's defects
    // Recombination only takes place within the cascade of an individual PKA.
    // Cascades are assumed to be non-overlapping in time and space.
    ///
    if (_trim_parameters.recombination && !_vacancy_buffer.empty())
    {
      // 1. build kd-tree for the vacancies
      const unsigned int max_leaf_size = 50; // slightly affects runtime
      auto point_list = PointListAdaptor<MyTRIMDefectBufferItem>(_vacancy_buffer.begin(), _vacancy_buffer.end());
      auto kd_tree = libmesh_make_unique<KDTreeType>(
          LIBMESH_DIM, point_list,
          nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));

      mooseAssert(kd_tree != nullptr, "KDTree was not properly initialized.");
      kd_tree->buildIndex();

      const Real r_rec = _trim_parameters.r_rec;
      nanoflann::SearchParams params;

      // 2. iterate over interstitials and recombine them if they are with r_rec of a vacancy
      std::vector<std::pair<std::size_t, Real>> ret_matches;
      for (auto & i: _interstitial_buffer)
      {
        ret_matches.clear();
        std::size_t n_result = kd_tree->radiusSearch(&(i.first(0)), r_rec, ret_matches, params);

        for (std::size_t j = 0; j < n_result; ++j)
        {
          auto & v = _vacancy_buffer[ret_matches[j].first];

          // only allow interstitial to go into vacancy of the same type
          if (v.second == i.second)
          {
            // mark vacancy-interstitial pair for deletion
            i.second = libMesh::invalid_uint;
            v.second = libMesh::invalid_uint;
            break;
          }
        }
      }
    }


    // add remaining defects to result
    for (auto & i: _interstitial_buffer) // the should happen above
      if (i.second != libMesh::invalid_uint)
        addDefectToResult(_rasterizer.periodicPoint(i.first), i.second, INTERSTITIAL);
    for (auto & v: _vacancy_buffer)
      if (v.second != libMesh::invalid_uint)
        addDefectToResult(_rasterizer.periodicPoint(v.first), v.second, VACANCY);
  }
}
