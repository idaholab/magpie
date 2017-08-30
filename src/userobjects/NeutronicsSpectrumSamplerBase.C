/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "NeutronicsSpectrumSamplerBase.h"
#include "MooseMesh.h"

#ifdef RATTLESNAKE_ENABLED
  #include "YakxsUtilities.h"
#endif

// C++ includes
#include <sstream>
#include <algorithm>
#include <limits>

template<>
InputParameters validParams<NeutronicsSpectrumSamplerBase>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredParam<std::vector<std::string> >("target_isotope_names", "The list of target isotope names e.g. U235.");
  params.addRequiredCoupledVar("number_densities", "Number densities for each isotope.");
  params.addRequiredParam<std::vector<Real> >("energy_group_boundaries", "The energy group boundaries in units of eV orderd "
                                                                         "from high to low [natural neutronics ordering].");
  params.addRequiredParam<unsigned int>("L", "The maximum order of Legendre terms for expanding the recoil XS in mu_lab.");
  params.addRequiredParam<std::vector<Point> >("points", "The points where you want to evaluate the variables");
  params.addClassDescription("Radiation Damage user object base class.\n Computes PDFs from neutronics calculations that are used to sample PKAs in BCMC simulations.");
  return params;
}

NeutronicsSpectrumSamplerBase::NeutronicsSpectrumSamplerBase(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _target_isotope_names(getParam<std::vector<std::string> >("target_isotope_names")),
    _energy_group_boundaries(getParam<std::vector<Real> >("energy_group_boundaries")),
    _I(_target_isotope_names.size()),
    _G(_energy_group_boundaries.size() - 1),
    _L(getParam<unsigned int>("L")),
    _points(getParam<std::vector<Point> >("points")),
    _npoints(_points.size()),
    _qp_is_cached(false)
{
  // check input dimensions
  if (coupledComponents("number_densities") != _I)
    mooseError("ZAID and number_densities must have the same length.");

  // get the number density variables
  _number_densities.resize(_I);
  for (unsigned int i = 0; i < _I; ++i)
    _number_densities[i] = & coupledValue("number_densities", i);

  _zaids.resize(_I);
  for (unsigned int i = 0; i < _I; ++i)
  {
#ifdef RATTLESNAKE_ENABLED
    unsigned int Z, A;
    // check if isotope names are valid, NOTE: error handling is delegated to Yakxs::Utilities
    YAKXS::Utility::getAZFromIsotopeName(_target_isotope_names[i], A, Z);
    // convert from name to ZAID
    _zaids[i] = YAKXS::Utility::stringToZaid(_target_isotope_names[i]);
#else
    _zaids[i] = localStringToZaid(_target_isotope_names[i]);
#endif
  }

  _owner.resize(_npoints);
  _qp_cache.resize(_npoints);

  // check energy group ordering
  if (_energy_group_boundaries.size() < 2)
    mooseError("At least 2 boundaries must be provided for energy_group_boundaries");
  for (unsigned int j = 1; j < _energy_group_boundaries.size(); ++j)
    if (_energy_group_boundaries[j] >= _energy_group_boundaries[j - 1])
      mooseError("energy_group_boundaries must be ordered from high to low");
}

void
NeutronicsSpectrumSamplerBase::execute()
{
  if (_local_elem_to_contained_points.find(_current_elem) != _local_elem_to_contained_points.end())
  {
    for (auto & j : _local_elem_to_contained_points[_current_elem])
    {
      _current_point = j;
      // NOTE: PKA is evaluated at a point that is not necessarily a qp but material
      // props only live on qps. This is a cheap and dirty solution evaluating the
      // PKA at the closest quadrature point.
      if (!_qp_is_cached)
      {
        Real min_dsq = (_q_point[0] - _points[_current_point]).norm_sq();
        unsigned int min_qp = 0;
        for (unsigned int qp = 1; qp < _q_point.size(); qp++)
        {
          Real dsq = (_q_point[qp] - _points[_current_point]).norm_sq();
          if (dsq < min_dsq)
          {
            min_dsq = dsq;
            min_qp = qp;
          }
        }
        _qp_cache[_current_point] = min_qp;
      }

      // Set current qp to cached min_qp
      _qp = _qp_cache[_current_point];

      // call back before computing radiation damage PDF for caching things that
      // don't change
      preComputeRadiationDamagePDF();

      // Evalute the radiation damage PDF at min_qp
      for (unsigned int i = 0; i < _I; ++i)
        for (unsigned int g = 0; g < _G; ++g)
          for (unsigned int p = 0; p < _nphi; ++p)
            for (unsigned int q = 0; q < _nmu; ++q)
              _sample_point_data[_current_point]({i, g, p, q}) = computeRadiationDamagePDF(i, g, p, q);
    }
  }
}

void
NeutronicsSpectrumSamplerBase::preComputeRadiationDamagePDF()
{
}

void
NeutronicsSpectrumSamplerBase::meshChanged()
{
  // make sure that _local_elem_to_contained_points is empty
  _local_elem_to_contained_points.clear();

  // Rebuild point_locator & find the right element for each point
  UniquePtr<PointLocatorBase> point_locator = _mesh.getPointLocator();
  for (unsigned int j = 0; j < _npoints; j++)
  {
    // make sure element is local
    _owner[j] = 0;
    const Elem * candidate_elem = (*point_locator)(_points[j]);
    if (candidate_elem && candidate_elem->processor_id() == processor_id())
    {
      _owner[j] = processor_id();
      if (_local_elem_to_contained_points.find(candidate_elem) != _local_elem_to_contained_points.end())
        _local_elem_to_contained_points[candidate_elem].push_back(j);
      else
        _local_elem_to_contained_points[candidate_elem] = {j};
    }
  }

  // make sure everyone knows who owns which point
  for (unsigned int j = 0; j < _npoints; j++)
    gatherMax(_owner[j]);

  _qp_is_cached = false;
}

void
NeutronicsSpectrumSamplerBase::initialSetup()
{
  // allocate PDF
  // NOTE: Needs to be delayed to initialize because _nSH is set in derived class
  for (unsigned int j = 0; j < _npoints; ++j)
    _sample_point_data.push_back(MultiIndex<Real>({_I, _G, _nphi, _nmu}));
}

void
NeutronicsSpectrumSamplerBase::initialize()
{
  meshChanged();

  for (unsigned j = 0; j < _npoints; ++j)
    for (auto entry : _sample_point_data[j])
      entry.second = 0;
}

void
NeutronicsSpectrumSamplerBase::finalize()
{
  if (_mesh.n_processors() > 1)
  {
    for (unsigned j = 0; j < _npoints; ++j)
    {
      std::vector<Real> flat_data(_I * _G * _nmu * _nphi);
      if (_owner[j] == processor_id())
        flat_data = _sample_point_data[j].getRawData();

      _communicator.broadcast(flat_data, _owner[j]);

      if (_owner[j] != processor_id())
        _sample_point_data[j] = MultiIndex<Real>({_I, _G, _nphi, _nmu}, flat_data);
    }
  }

  // set _qp_is_cached flag to true. Actually we only need to do this if
  // it was false but this way we can save an if statement
  _qp_is_cached = true;
}

void
NeutronicsSpectrumSamplerBase::threadJoin(const UserObject & y)
{
  const NeutronicsSpectrumSamplerBase & uo = static_cast<const NeutronicsSpectrumSamplerBase &>(y);
  for (unsigned j = 0; j < _npoints; ++j)
    for (unsigned int i = 0; i < _I; ++i)
      for (unsigned int g = 0; g < _G; ++g)
        for (unsigned int p = 0; p < _nphi; ++p)
          for (unsigned int q = 0; q < _nmu; ++q)
            _sample_point_data[j]({i, g, p, q}) += uo._sample_point_data[j]({i, g, p, q});
}

MultiIndex<Real>
NeutronicsSpectrumSamplerBase::getPDF(unsigned int point_id) const
{
  return _sample_point_data[point_id];
}

std::vector<unsigned int>
NeutronicsSpectrumSamplerBase::getZAIDs() const
{
  return _zaids;
}

std::vector<Real>
NeutronicsSpectrumSamplerBase::getEnergies() const
{
  return _energy_group_boundaries;
}

bool
NeutronicsSpectrumSamplerBase::hasIsotope(std::string target_isotope) const
{
  auto it = std::find(_target_isotope_names.begin(), _target_isotope_names.end(), target_isotope);
  return it != _target_isotope_names.end();
}

unsigned int
NeutronicsSpectrumSamplerBase::localStringToZaid(std::string s) const
{
  if (s == "U235")
    return 922350;
  else if (s == "U238")
    return 922380;
  else if (s == "Pu238")
    return 942380;
  else if (s == "Pu239")
    return 942390;
  else if (s == "Pu240")
    return 942400;
  mooseError("Isotope name ", s, " cannot be converted with the localStringToZaid conversion method.");
  return 0;
}
