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
PointListAdaptor<ThreadedRecoilLoopBase::MyTRIMDefectBufferItem>::getPoint(
    const ThreadedRecoilLoopBase::MyTRIMDefectBufferItem & item) const
{
  return item.point;
}

ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const MyTRIMRasterizer & rasterizer,
                                               const MooseMesh & mesh)
  : _rasterizer(rasterizer),
    _trim_parameters(_rasterizer.getTrimParameters()),
    _nvars(_trim_parameters.nVars()),
    _mesh(mesh),
    _dim(_mesh.dimension())
{
  _simconf.setLengthScale(_trim_parameters.length_scale);
}

// Splitting Constructor
ThreadedRecoilLoopBase::ThreadedRecoilLoopBase(const ThreadedRecoilLoopBase & x,
                                               Threads::split /*split*/)
  : _rasterizer(x._rasterizer),
    _trim_parameters(x._trim_parameters),
    _nvars(x._nvars),
    _mesh(x._mesh),
    _dim(x._dim)
{
  _simconf.setLengthScale(_trim_parameters.length_scale);
}

void
ThreadedRecoilLoopBase::operator()(const PKARange & pka_list)
{
  // fetch a point locator
  _pl = _mesh.getPointLocator();

  // permit querying points that are potentially outside the mesh
  _pl->enable_out_of_mesh_mode();

  // create a new sample class to bridge the MOOSE mesh and the MyTRIM domain
  MooseMyTRIMSample sample(_rasterizer, _mesh, &_simconf);

  // create a FIFO for recoils
  std::queue<MyTRIM_NS::IonBase *> recoils;

  // create a list potentially used for energy deposition
  std::list<std::pair<Point, Real>> edep_list;

  // build the requested TRIM module
  std::unique_ptr<MooseMyTRIMCore> TRIM;
  switch (_trim_parameters.trim_module)
  {
    // basic module with interstitial and vacancy generation
    case MyTRIMRasterizer::MYTRIM_CORE:
      TRIM.reset(new MooseMyTRIMCore(&_simconf, &sample));
      break;

    // record energy deposited to the lattice
    case MyTRIMRasterizer::MYTRIM_ENERGY_DEPOSITION:
      TRIM.reset(new MooseMyTRIMEnergyDeposition(&_simconf, &sample, edep_list));
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
    // seed the RNG with the seed assigned to the primary knock on atom (for parallel
    // reproducibility)
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

      // deal with what this recoil left behind
      if (recoil->_state == MyTRIM_NS::IonBase::VACANCY)
        _vacancy_buffer.push_back(MyTRIMDefectBufferItem(recoil->_pos, recoil->_tag));
      else if (recoil->_state == MyTRIM_NS::IonBase::REPLACEMENT)
        addDefectToResult(
            _rasterizer.periodicPoint(recoil->_pos), recoil->_tag, 1, REPLACEMENT_OUT);
      // ignore like-for-like substitutions

      // remove zero energy recoils
      if (recoil->_E > 0.0)
      {
        // full recoil or analytical approximation
        if (recoil->_E < _trim_parameters.analytical_cutoff &&
            _rasterizer.isTrackedSpecies(recoil->_Z, recoil->_m))
        {
          mooseDoOnce(mooseWarning("Skipping detailed cascade calculation below cutoff energy."));
          const auto pp = _rasterizer.periodicPoint(recoil->_pos);

#ifdef GSL_ENABLED
          const auto elem = (*_pl)(pp);

          // the actual number density is immaterial for NRT, only number fractions are important
          // so we go ahead and first get it from the rasterizer and then normalize it to sum to 1
          auto number_fractions = _rasterizer.material(elem);

          Real total_number_density = 0;
          for (unsigned int j = 0; j < number_fractions.size(); ++j)
            total_number_density += number_fractions[j];

          for (unsigned int j = 0; j < number_fractions.size(); ++j)
            number_fractions[j] /= total_number_density;

          unsigned int index;
          Real distance;

          findBestNRTMatch(number_fractions, index, distance);
          if (distance > _trim_parameters.max_nrt_distance)
            index = addNRTEntry(number_fractions);

          // look up entry for recoil Z & m
          unsigned int species_index = _pa_nrt[index]->findSpeciesIndex(recoil->_Z, recoil->_m);

          // interpolate the replacement counts for all considered target species
          // vacancy energy is the energy of the cascade that goes into creating vacancies
          // vacancy_energy = sum_{targets: j} Ebind_j * n_{projectile, j}
          Real vacancy_energy = 0;
          if (species_index != libMesh::invalid_uint)
            for (unsigned int target_var = 0; target_var < _nvars; ++target_var)
            {
              unsigned int target_species_index = _pa_nrt[index]->findSpeciesIndex(
                  _trim_parameters.element_prototypes[target_var]._Z,
                  _trim_parameters.element_prototypes[target_var]._m);
              // using linear perturbation theory estimate of g_ij(number_fractions)
              Real w = _pa_nrt[index]->linearInterpolation(
                  species_index, target_species_index, recoil->_E);
              for (unsigned int l = 0; l < _pa_derivative_nrt[index]->nSpecies(); ++l)
                w += _pa_derivative_nrt[index]->linearInterpolation(
                         species_index, target_species_index, l, recoil->_E) *
                     (number_fractions[l] - _pa_nrt[index]->numberFraction(l));

              // increment energy, vacancy and interstitial buffers
              vacancy_energy += w * _trim_parameters.element_prototypes[target_var]._Elbind;
              _vacancy_buffer.emplace_back(MyTRIMDefectBufferItem(recoil->_pos, target_var, w));
              _interstitial_buffer.emplace_back(
                  MyTRIMDefectBufferItem(recoil->_pos, target_var, w));
            }
          else
            mooseError("NRT treatment did not find recoil with ", recoil->_Z, " ", recoil->_m);

          // add remaining recoil energy
          if (_trim_parameters.trim_module == MyTRIMRasterizer::MYTRIM_ENERGY_DEPOSITION)
            addEnergyToResult(pp, recoil->_E - vacancy_energy);
#else
          mooseDoOnce(
              mooseWarning("GSL is not enabled and cutoff energy != 0. Polyatomic NRT requires "
                           "GSL, current settings do not account for skipped ions."));
          // add remaining recoil energy
          if (_trim_parameters.trim_module == MyTRIMRasterizer::MYTRIM_ENERGY_DEPOSITION)
            addEnergyToResult(pp, recoil->_E);
#endif
        }
        else
        {
          // follow this ion's trajectory and store recoils
          TRIM->trim(recoil, recoils);

          // are we tracking atoms of this type?
          if (recoil->_tag >= 0)
          {
            if (recoil->_state == MyTRIM_NS::IonBase::INTERSTITIAL)
              _interstitial_buffer.emplace_back(MyTRIMDefectBufferItem(recoil->_pos, recoil->_tag));
            else if (recoil->_state == MyTRIM_NS::IonBase::REPLACEMENT)
              addDefectToResult(
                  _rasterizer.periodicPoint(recoil->_pos), recoil->_tag, 1, REPLACEMENT_IN);
            // BUG: untracked PKAs in replacement collisions will cause a REPLACEMENT_OUT without a
            // REPLACEMENT_IN!
          }

          // store energy deposition
          for (auto & edep : edep_list)
            addEnergyToResult(_rasterizer.periodicPoint(edep.first), edep.second);
          edep_list.clear();
        }
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
      auto point_list =
          PointListAdaptor<MyTRIMDefectBufferItem>(_vacancy_buffer.begin(), _vacancy_buffer.end());
      auto kd_tree = libmesh_make_unique<KDTreeType>(
          LIBMESH_DIM, point_list, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));

      mooseAssert(kd_tree != nullptr, "KDTree was not properly initialized.");
      kd_tree->buildIndex();

      const Real r_rec2 = _trim_parameters.r_rec * _trim_parameters.r_rec;
      nanoflann::SearchParams params;

      // 2. iterate over interstitials and recombine them if they are with r_rec of a vacancy
      std::vector<std::pair<std::size_t, Real>> ret_matches;
      for (auto & i : _interstitial_buffer)
      {
        ret_matches.clear();
        std::size_t n_result = kd_tree->radiusSearch(&(i.point(0)), r_rec2, ret_matches, params);

        for (std::size_t j = 0; j < n_result; ++j)
        {
          auto & v = _vacancy_buffer[ret_matches[j].first];

          // only allow interstitial to go into vacancy of the same type
          if (v.variable_id == i.variable_id)
          {
            // mark vacancy-interstitial pair for deletion
            i.variable_id = libMesh::invalid_uint;
            v.variable_id = libMesh::invalid_uint;
            break;
          }
        }
      }
    }

    // add remaining defects to result
    for (auto & i : _interstitial_buffer)
      if (i.variable_id != libMesh::invalid_uint)
        addDefectToResult(
            _rasterizer.periodicPoint(i.point), i.variable_id, i.weight, INTERSTITIAL);
    for (auto & v : _vacancy_buffer)
      if (v.variable_id != libMesh::invalid_uint)
        addDefectToResult(_rasterizer.periodicPoint(v.point), v.variable_id, v.weight, VACANCY);
  }
}

#ifdef GSL_ENABLED
unsigned int
ThreadedRecoilLoopBase::addNRTEntry(const std::vector<Real> & number_fractions)
{
  if (number_fractions.size() != _nvars)
    mooseError("number_fractions has wrong size");

  std::vector<MyTRIM_NS::Element> poly_mat;
  for (unsigned int j = 0; j < _nvars; ++j)
  {
    MyTRIM_NS::Element element;
    element._Z = _trim_parameters.element_prototypes[j]._Z;
    element._m = _trim_parameters.element_prototypes[j]._m;
    element._t = number_fractions[j];
    element._Edisp = _trim_parameters.element_prototypes[j]._Edisp;
    element._Elbind = _trim_parameters.element_prototypes[j]._Elbind;
    poly_mat.push_back(element);
  }

  /**
   * NOTE: using the net displacement rate here defined as
   *       "[...] atoms displaced and not recaptured in subsequent replacement collisions"
   * TODO: this should be checked for accuracy, also net displacement rate is not as well
   *       verified as total displacement rate
   */
  _pa_nrt.push_back(libmesh_make_unique<PolyatomicDisplacementFunction>(poly_mat, NET));
  _pa_derivative_nrt.push_back(libmesh_make_unique<PolyatomicDisplacementDerivativeFunction>(
      poly_mat, NET_DERIVATIVE, _pa_nrt.back().get()));

  // integrate P&K's equation to the analytical cutoff
  Real energy = _pa_nrt.back()->minEnergy();

  // some baby steps [10 eV] to get initial evolution right
  Real upper_linear_range = std::max(2 * energy, 100.0);
  unsigned int n_initial_steps = std::ceil((upper_linear_range - energy) / 10);
  for (unsigned int j = 0; j < n_initial_steps; ++j)
  {
    energy += 10;
    _pa_nrt.back()->advanceDisplacements(energy);
  }

  // constant in log(E) stepping
  for (;;)
  {
    energy *= _trim_parameters.nrt_log_energy_spacing;
    if (energy > _trim_parameters.analytical_cutoff)
    {
      _pa_nrt.back()->advanceDisplacements(_trim_parameters.analytical_cutoff);
      break;
    }
    _pa_nrt.back()->advanceDisplacements(energy);
  }

  for (unsigned int n = 1; n < _pa_nrt.back()->nEnergySteps(); ++n)
    _pa_derivative_nrt.back()->advanceDisplacements(_pa_nrt.back()->energyPoint(n));

  return _pa_nrt.size() - 1;
}

void
ThreadedRecoilLoopBase::findBestNRTMatch(const std::vector<Real> & number_fractions,
                                         unsigned int & index,
                                         Real & distance) const
{
  index = 0;
  distance = std::numeric_limits<Real>::max();

  if (_pa_nrt.size() == 0)
    return;

  for (unsigned int i = 0; i < _pa_nrt.size(); ++i)
  {
    auto & nrt_candidate = _pa_nrt[i];

    mooseAssert(number_fractions.size() == nrt_candidate->nSpecies(),
                "Number densities have different sizes");

    // find the infinity norm difference
    Real inf_norm_difference = 0;
    for (unsigned int j = 0; j < number_fractions.size(); ++j)
      if (std::abs(number_fractions[j] - nrt_candidate->numberFraction(j)) > inf_norm_difference)
        inf_norm_difference = std::abs(number_fractions[j] - nrt_candidate->numberFraction(j));

    if (distance > inf_norm_difference)
    {
      distance = inf_norm_difference;
      index = i;
    }
  }
}
#endif
