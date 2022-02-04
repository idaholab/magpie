/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMRasterizer.h"
#include "PKAGeneratorBase.h"
#include "MooseRandom.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"

#include <type_traits>

// custom data load and data store methods for a struct with an std::vector member
template <>
inline void
dataStore(std::ostream & stream, MyTRIMRasterizer::AveragedData & ad, void * context)
{
  dataStore(stream, ad._elements, context);
  dataStore(stream, ad._site_volume, context);
}
template <>
inline void
dataLoad(std::istream & stream, MyTRIMRasterizer::AveragedData & ad, void * context)
{
  dataLoad(stream, ad._elements, context);
  dataLoad(stream, ad._site_volume, context);
}

// custom data load and data store methods for a class with virtual members (vtable pointer must not
// be (un)serialized)
template <>
inline void
dataStore(std::ostream & stream, MyTRIM_NS::IonBase & pl, void * context)
{
  dataStore(stream, pl._Z, context);
  dataStore(stream, pl._m, context);
  dataStore(stream, pl._E, context);
  dataStore(stream, pl._dir, context);
  dataStore(stream, pl._pos, context);
  dataStore(stream, pl._seed, context);
  dataStore(stream, pl._gen, context);
  dataStore(stream, pl._id, context);
  dataStore(stream, pl._tag, context);
  dataStore(stream, pl._Ef, context);
  dataStore(stream, pl._state, context);
}
template <>
inline void
dataLoad(std::istream & stream, MyTRIM_NS::IonBase & pl, void * context)
{
  dataLoad(stream, pl._Z, context);
  dataLoad(stream, pl._m, context);
  dataLoad(stream, pl._E, context);
  dataLoad(stream, pl._dir, context);
  dataLoad(stream, pl._pos, context);
  dataLoad(stream, pl._seed, context);
  dataLoad(stream, pl._gen, context);
  dataLoad(stream, pl._id, context);
  dataLoad(stream, pl._tag, context);
  dataLoad(stream, pl._Ef, context);
  dataLoad(stream, pl._state, context);
}

registerMooseObject("MagpieApp", MyTRIMRasterizer);

InputParameters
MyTRIMRasterizer::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription("Gather the element distribution of the simulation domain for a TRIM "
                             "binary collision Monte Carlo simulation");

  MooseEnum var_physical_meaning("STOICHIOMETRY NUMBER_DENSITY", "STOICHIOMETRY");
  params.addParam<MooseEnum>("var_physical_meaning",
                             var_physical_meaning,
                             "The physical meaning of the rasterizer variables.");

  params.addRequiredCoupledVar("var", "Variables to rasterize");
  params.addCoupledVar("periodic_var",
                       "Optional variables that determines the periodicity. If not supplied the "
                       "first argument of 'var' will be used.");
  params.addRequiredParam<std::vector<Real>>("M", "Element mass in amu");
  params.addRequiredParam<std::vector<Real>>("Z", "Nuclear charge in e");
  params.addParam<std::vector<Real>>("Mtol",
                                     "Tolerance on mass number for tagging PKAs with var id.");
  params.addParam<std::vector<Real>>("Ebind", "Binding energy in eV");
  params.addParam<std::vector<Real>>("Edisp", "Displacement threshold in eV");
  params.addParam<MaterialPropertyName>(
      "site_volume", "Lattice site volume in nm^3 (regardless of the chosen mesh units)");
  params.addRequiredParam<std::vector<UserObjectName>>("pka_generator",
                                                       "List of PKA generating user objects");
  ExecFlagEnum setup_options(MooseUtils::getDefaultExecFlagEnum());

  // we run this object once a timestep
  setup_options = EXEC_TIMESTEP_BEGIN;
  params.set<ExecFlagEnum>("execute_on") = setup_options;

  // which TRIM Module to run for optional capabilities like energy deposition
  MooseEnum trim_module_options("CORE=0 ENERGY_DEPOSITION=1", "CORE");
  params.addParam<MooseEnum>("trim_module",
                             trim_module_options,
                             "TRIM Module to run for optional capabilities like energy deposition");

  // which units of length to use
  MooseEnum length_unit_options("ANGSTROM=0 NANOMETER MICROMETER", "ANGSTROM");
  params.addParam<MooseEnum>("length_unit",
                             length_unit_options,
                             "Length units of the MOOSE mesh. MyTRIM contains pretabulated "
                             "crossection data with units so this option must be set correctly to "
                             "obtain physical results.");

  // Advanced options
  params.addParam<unsigned int>("interval", 1, "The time step interval at which TRIM BCMC is run");
  params.addParam<Real>("analytical_energy_cutoff",
                        0.0,
                        "Energy cutoff in eV below which recoils are not followed explicitly but "
                        "effects are calculated analytically.");
  params.addParam<Real>("recoil_rate_scaling",
                        1.0,
                        "A factor to scale computed reaction rates in the the PKAGenerator "
                        "objects. This is useful to avoid extremely large PKA lists.");
  params.addParam<unsigned int>(
      "max_pka_count", "Desired number of PKAs to be run during each invocation of mytrim");
  params.addRangeCheckedParam<Real>("max_nrt_difference",
                                    0.2,
                                    "max_nrt_difference > 0 & max_nrt_difference < 1",
                                    "The largest max-norm difference between number fractions for "
                                    "reusing existing polyatomic NRT.");
  params.addParamNamesToGroup("interval analytical_energy_cutoff max_pka_count", "Advanced");

  params.addParam<Real>("r_rec",
                        "Recombination radius in Angstrom. Frenkel pairs with a separation "
                        "distance lower than this will be removed from the cascade");
  params.addParamNamesToGroup("r_rec", "Recombination");

  return params;
}

MyTRIMRasterizer::MyTRIMRasterizer(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _nvars(coupledComponents("var")),
    _dim(_mesh.dimension()),
    _var(_nvars),
    _site_volume_prop(nullptr),
    _pka_generator_names(getParam<std::vector<UserObjectName>>("pka_generator")),
    _pka_generators(),
    _periodic(isCoupled("periodic_var") ? coupled("periodic_var", 0) : coupled("var", 0)),
    _accumulated_time(0.0),
    _accumulated_time_old(0.0),
    _interval(getParam<unsigned int>("interval")),
    _perf_finalize(registerTimedSection("finalize", 2))
{
  if (_nvars == 0)
    mooseError("Must couple variables to MyTRIMRasterizer.");

  // fill in trim parameters
  _trim_parameters.element_prototypes.resize(_nvars);
  _trim_parameters.analytical_cutoff = getParam<Real>("analytical_energy_cutoff");
  _trim_parameters.recoil_rate_scaling = getParam<Real>("recoil_rate_scaling");
  _trim_parameters.max_nrt_distance = getParam<Real>("max_nrt_difference");
  _trim_parameters.nrt_log_energy_spacing = 1.1;
  _trim_parameters.trim_module = getParam<MooseEnum>("trim_module").getEnum<TRIMModuleEnum>();
  if (isParamValid("max_pka_count"))
    _trim_parameters.desired_npka = getParam<unsigned int>("max_pka_count");
  else
    _trim_parameters.desired_npka = 0;

  auto trim_M = getParam<std::vector<Real>>("M");
  auto trim_Z = getParam<std::vector<Real>>("Z");

  std::vector<Real> mtol;
  if (isParamValid("Mtol"))
    mtol = getParam<std::vector<Real>>("Mtol");
  else
    mtol.assign(_nvars, 0.5);

  if (trim_M.size() != _nvars)
    mooseError("Parameter 'M' must have as many components as coupled variables.");
  if (trim_Z.size() != _nvars)
    mooseError("Parameter 'Z' must have as many components as coupled variables.");
  if (mtol.size() != _nvars)
    mooseError("Parameter 'mtol' must have as many components as coupled variables.");

  // error check masses and charges
  for (unsigned int i = 0; i < _nvars; ++i)
  {
    if (trim_Z[i] > trim_M[i])
      mooseError("Value of Z is larger than value of M for entry ", i);
    if (trim_Z[i] >= _pka_parameters._index_Z.size())
      mooseError("Value of Z is too large. Maximum Z supported is ",
                 _pka_parameters._index_Z.size() - 1,
                 " but one element has Z=",
                 i);
  }

  for (unsigned int i = 0; i < _nvars; ++i)
  {
    _var[i] = &coupledValue("var", i);
    _trim_parameters.element_prototypes[i]._m = trim_M[i];
    _trim_parameters.element_prototypes[i]._Z = trim_Z[i];
  }

  auto trim_Ebind = getParam<std::vector<Real>>("Ebind");
  if (trim_Ebind.size() == _nvars)
    for (unsigned int i = 0; i < _nvars; ++i)
      _trim_parameters.element_prototypes[i]._Elbind = trim_Ebind[i];
  else if (trim_Ebind.empty())
    for (unsigned int i = 0; i < _nvars; ++i)
      _trim_parameters.element_prototypes[i]._Elbind = 3.0;
  else
    mooseError("Parameter 'Ebind' must have as many components as coupled variables (or left empty "
               "for a default of 3eV).");

  auto trim_Edisp = getParam<std::vector<Real>>("Edisp");
  if (trim_Edisp.size() == _nvars)
    for (unsigned int i = 0; i < _nvars; ++i)
      _trim_parameters.element_prototypes[i]._Edisp = trim_Edisp[i];
  else if (trim_Edisp.empty())
    for (unsigned int i = 0; i < _nvars; ++i)
      _trim_parameters.element_prototypes[i]._Edisp = 25.0;
  else
    mooseError("Parameter 'Edisp' must have as many components as coupled variables (or left empty "
               "for a default of 25eV).");

  if (isParamValid("r_rec"))
  {
    _trim_parameters.recombination = true;
    _trim_parameters.r_rec = getParam<Real>("r_rec");
  }
  else
    _trim_parameters.recombination = false;

  // fetch PKA Generators
  for (auto && name : _pka_generator_names)
    _pka_generators.push_back(&getUserObjectByName<PKAGeneratorBase>(name));

  // set up data for sample periodicity
  for (unsigned int i = 0; i < _dim; ++i)
  {
    _pbc[i] = _mesh.isRegularOrthogonal() && _mesh.isTranslatedPeriodic(_periodic, i);

    if (_pbc[i])
    {
      _min_dim(i) = _mesh.getMinInDimension(i);
      _max_dim(i) = _mesh.getMaxInDimension(i);
    }
  }

  // determine length scale factor for TRIM
  switch (getParam<MooseEnum>("length_unit").getEnum<Unit>())
  {
    case ANGSTROM:
      _trim_parameters.length_scale = 1.0;
      break;

    case NANOMETER:
      _trim_parameters.length_scale = 10.0;
      break;

    case MICROMETER:
      _trim_parameters.length_scale = 10000.0;
      break;

    default:
      mooseError("Unknown length unit.");
  }

  if (getParam<MooseEnum>("var_physical_meaning") == "STOICHIOMETRY")
  {
    if (!isParamValid("site_volume"))
      mooseError("Rasterizer variables are stoiciometric contents, site_volume must be provided.");
    _site_volume_prop = &getMaterialProperty<Real>("site_volume");
  }
  else
    _site_volume_conversion = Utility::pow<3>(_trim_parameters.length_scale) * 1e-3;

  // setup invariant PKA generation parameters
  for (auto & nZ : _pka_parameters._index_Z)
    nZ = std::make_pair(0, 0);
  _pka_parameters._mass_charge_tuple.resize(_nvars);
  _pka_parameters._recoil_rate_scaling = _trim_parameters.recoil_rate_scaling;
  for (unsigned int i = 0; i < _nvars; ++i)
  {
    const auto Z = _trim_parameters.element_prototypes[i]._Z;

    // insert (mass, charge) pair
    _pka_parameters._mass_charge_tuple[i] =
        std::make_tuple(_trim_parameters.element_prototypes[i]._m, Z, mtol[i]);

    // increase the count of elements with the same Z
    auto & index_Z = _pka_parameters._index_Z[Z];
    index_Z.first++;

    // only set this to the first index (important for ionTag())
    if (index_Z.first == 1)
      index_Z.second = i;
  }
}

bool
MyTRIMRasterizer::executeThisTimestep() const
{
  return (_fe_problem.timeStep() - 1) % _interval == 0;
}

void
MyTRIMRasterizer::initialize()
{
  _execute_this_timestep = executeThisTimestep();

  // We roll back the accumulated time time if the preceeding timestep did
  // not converge
  if (!_fe_problem.converged())
    _accumulated_time = _accumulated_time_old;

  if (_execute_this_timestep)
  {
    _trim_parameters.last_executed_dt = _fe_problem.dt();
    _material_map.clear();
    _pka_list.clear();
  }

  /// setup global PKA parameters for the current timestep
  _pka_parameters._dt = _accumulated_time + _fe_problem.dt();
}

void
MyTRIMRasterizer::execute()
{
  // bail out early if not executing this timestep
  if (!_execute_this_timestep)
    return;

  // average element concentrations

  AveragedData average(_nvars);
  Real vol = 0.0;

  // average material data over elements
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
  {
    const Real qpvol = _JxW[qp] * _coord[qp];
    vol += qpvol;

    // average compositions on the element
    for (unsigned int i = 0; i < _nvars; ++i)
      average._elements[i] += qpvol * (*_var[i])[qp];

    // average site volume property
    if (_site_volume_prop)
      average._site_volume += qpvol * (*_site_volume_prop)[qp];
    else
      average._site_volume += qpvol * _site_volume_conversion;
  }

  // divide by total element volume
  if (vol > 0.0)
  {
    for (unsigned int i = 0; i < _nvars; ++i)
      average._elements[i] /= vol;

    average._site_volume /= vol;
  }

  // store in map
  _material_map[_current_elem->id()] = average;

  // update corrent element volume
  _pka_parameters._volume = vol;

  // add PKAs for current element
  for (auto && gen : _pka_generators)
    gen->appendPKAs(_pka_list, _pka_parameters, average);
}

void
MyTRIMRasterizer::threadJoin(const UserObject & y)
{
  // if the map needs to be updated we merge the maps from all threads
  if (_execute_this_timestep)
  {
    const MyTRIMRasterizer & uo = static_cast<const MyTRIMRasterizer &>(y);
    _material_map.insert(uo._material_map.begin(), uo._material_map.end());
    _pka_list.insert(_pka_list.end(), uo._pka_list.begin(), uo._pka_list.end());
  }
}

void
MyTRIMRasterizer::finalize()
{
  TIME_SECTION(_perf_finalize);

  // save the accumulated time so that we can properly roll back if the step does not converge
  _accumulated_time_old = _accumulated_time;

  // bail out early if not executing this timestep
  if (!_execute_this_timestep)
  {
    // no BCMC done, so wee accumulate this step's time to be taken care of by a later BCMC run
    _accumulated_time += _fe_problem.dt();
    return;
  }

  // BCMC was run for the accumulated time - the debt is paid
  _accumulated_time = 0.0;

  // for single processor runs we do not need to do anything here
  if (_app.n_processors() > 1)
  {
    // create send buffer
    std::string send_buffer;

    // create byte buffers for the streams received from all processors
    std::vector<std::string> recv_buffers;

    // pack the complex datastructures into the string stream
    serialize(send_buffer);

    // broadcast serialized data to and receive from all processors
    _communicator.allgather(send_buffer, recv_buffers);

    // unpack the received data and merge it into the local data structures
    deserialize(recv_buffers);
  }

  // we will assign random seeds on proc 0 & reject on proc 0. To guarantee reproducibility
  // we need to sort the PKA list
  if (processor_id() == 0)
  {
    std::vector<unsigned int> pka_seeds(_pka_list.size());

    // only do this on proc 0 thread 0
    for (auto i = beginIndex(_pka_list); i < _pka_list.size(); ++i)
      pka_seeds[i] = MooseRandom::randl();

    // sort PKA list only on processor 0 & assign random number seeds
    std::sort(_pka_list.begin(),
              _pka_list.end(),
              [](MyTRIM_NS::IonBase a, MyTRIM_NS::IonBase b)
              {
                return (a._pos < b._pos) || (a._pos == b._pos && a._m < b._m) ||
                       (a._pos == b._pos && a._m == b._m && a._E < b._E) ||
                       (a._pos == b._pos && a._m == b._m && a._E == b._E && a._Z < b._Z);
              });

    // store seeds in tag values
    for (auto i = beginIndex(_pka_list); i < _pka_list.size(); ++i)
      _pka_list[i]._seed = pka_seeds[i];
  }

  // rejection is performed in processor 0 only
  _trim_parameters.original_npka = _pka_list.size();
  if (_trim_parameters.desired_npka == 0 ||
      _trim_parameters.desired_npka > _trim_parameters.original_npka)
  {
    _trim_parameters.scaled_npka = _trim_parameters.original_npka;
    _trim_parameters.result_scaling_factor = 1.0 / _trim_parameters.recoil_rate_scaling;
  }
  else
  {
    if (processor_id() == 0)
    {
      Real acceptance_probability =
          Real(_trim_parameters.desired_npka) / Real(_trim_parameters.original_npka);

      // most straight-forward but probably inefficient implementation of rejection
      std::vector<MyTRIM_NS::IonBase> old_pka_list = _pka_list;
      _pka_list.resize(0);
      for (auto & p : old_pka_list)
        if (MooseRandom::rand() < acceptance_probability)
          _pka_list.push_back(p);

      // save the size of the PKA list after rejection & the scaling factor
      _trim_parameters.scaled_npka = _pka_list.size();
      _trim_parameters.result_scaling_factor = Real(_trim_parameters.original_npka) /
                                               Real(_trim_parameters.scaled_npka) /
                                               _trim_parameters.recoil_rate_scaling;
    }

    // need to broadcast the size of the PKA list after rejection & result scaling factor
    _communicator.broadcast(_trim_parameters.scaled_npka);
    _communicator.broadcast(_trim_parameters.result_scaling_factor);
  }

  // communicate the PKA list if n_proc > 1
  if (_app.n_processors() > 1)
  {
    std::string pka_list_buffer;
    if (processor_id() == 0)
    {
      // pack the local _pka_list into a string buffer
      std::ostringstream oss;
      dataStore(oss, _pka_list, this);
      pka_list_buffer.assign(oss.str());
    }

    // communicate the pka list
    _communicator.broadcast(pka_list_buffer);
    if (processor_id() != 0)
    {
      _pka_list.resize(0);
      std::istringstream iss(pka_list_buffer);
      dataLoad(iss, _pka_list, this);
    }
  }

  // prune PKA list
  if (_app.n_processors() > 1)
  {
    // split PKAs into per-processor ranges
    std::vector<unsigned int> interval(_app.n_processors() + 1, 0);
    for (unsigned int i = 0; i < _app.n_processors(); ++i)
      interval[i + 1] = (_pka_list.size() - interval[i]) / (_app.n_processors() - i) + interval[i];

    auto begin = interval[processor_id()];
    auto end = interval[processor_id() + 1];
    std::vector<MyTRIM_NS::IonBase> own_pka_list(end - begin);

    for (auto i = begin; i < end; ++i)
      own_pka_list[i - begin] = _pka_list[i];

    _pka_list.swap(own_pka_list);
  }
}

const std::vector<Real> &
MyTRIMRasterizer::material(const Elem * elem) const
{
  auto i = _material_map.find(elem->id());

  // there should be data for every element in the mesh
  if (i == _material_map.end())
    mooseError("Element not found in material map.");

  return i->second._elements;
}

Real
MyTRIMRasterizer::siteVolume(const Elem * elem) const
{
  auto i = _material_map.find(elem->id());

  // there should be data for every element in the mesh
  if (i == _material_map.end())
    mooseError("Element not found in material map.");

  return i->second._site_volume;
}

Point
MyTRIMRasterizer::periodicPoint(const Point & pos) const
{
  // point to sample the material at
  Point p(pos(0), pos(1), _dim == 2 ? 0.0 : pos(2));

  // apply periodic boundary conditions
  for (unsigned int i = 0; i < _dim; ++i)
    if (_pbc[i])
    {
      const Real width = _max_dim(i) - _min_dim(i);
      p(i) -= std::floor((p(i) - _min_dim(i)) / width) * width;
    }

  return p;
}

bool
MyTRIMRasterizer::isTrackedSpecies(unsigned int atomic_number, Real mass_number) const
{
  mooseAssert(_pka_generators[0], "PKA generator is not set.");
  return _pka_generators[0]->ionTag(_pka_parameters, atomic_number, mass_number) != -1;
}

void
MyTRIMRasterizer::serialize(std::string & serialized_buffer)
{
  // stream for serializing the _material_map and _pka_list data structure to a byte stream
  std::ostringstream oss;
  dataStore(oss, _material_map, this);
  dataStore(oss, _pka_list, this);

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
}

void
MyTRIMRasterizer::deserialize(std::vector<std::string> & serialized_buffers)
{
  mooseAssert(serialized_buffers.size() == _app.n_processors(),
              "Unexpected size of serialized_buffers: " << serialized_buffers.size());

  // The input string stream used for deserialization
  std::istringstream iss;

  // Loop over all datastructures for all procerssors to perfrom the gather operation
  for (unsigned int rank = 0; rank < serialized_buffers.size(); ++rank)
  {
    // skip the current processor (its data is already in the structures)
    if (rank == processor_id())
      continue;

    // populate the stream with a new buffer and reset stream state
    iss.str(serialized_buffers[rank]);
    iss.clear();

    // Load the communicated data into temporary structures
    MaterialMap other_material_map;
    dataLoad(iss, other_material_map, this);
    std::vector<MyTRIM_NS::IonBase> other_pka_list;
    dataLoad(iss, other_pka_list, this);

    // merge the data in with the current processor's data
    _material_map.insert(other_material_map.begin(), other_material_map.end());

    // merging the PKA lists
    _pka_list.insert(_pka_list.begin(), other_pka_list.begin(), other_pka_list.end());
  }
}
