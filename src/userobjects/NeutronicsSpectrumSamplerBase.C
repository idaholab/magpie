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
  params.addRequiredParam<std::vector<Real> >("energy_group_boundaries", "The energy group boundaries ommitting E = 0.0. Units are MeV.");
  params.addRequiredParam<unsigned int>("L", "The order up to which angular moments of the PKA distribution are computed.");
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
    _zaids[i] = stringToZaid(_target_isotope_names[i]);
#else
    _zaids[i] = localStringToZaid(_target_isotope_names[i]);
#endif
  }

  _point_element.resize(_npoints);
  _qp_cache.resize(_npoints);
}

void
NeutronicsSpectrumSamplerBase::execute()
{
  if (std::find(_point_element.begin(), _point_element.end(), _current_elem->id()) != _point_element.end())
  {
    for (_current_point = 0; _current_point < _npoints; _current_point++)
    {
      if (_point_element[_current_point] != _current_elem->id()) continue;

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
          for (unsigned int p = 0; p < _nSH; ++p)
            _sample_point_data[_current_point]({i, g, p}) = computeRadiationDamagePDF(i, g, p);
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
  // Rebuild point_locator & find the right element for each point
  UniquePtr<PointLocatorBase> point_locator = PointLocatorBase::build(TREE_ELEMENTS, _mesh.getMesh());
  for (unsigned int j = 0; j < _npoints; j++)
    _point_element[j] = (*point_locator)(_points[j])->id();
  _qp_is_cached = false;
}

void
NeutronicsSpectrumSamplerBase::initialSetup()
{
  // allocate PDF
  // NOTE: Needs to be delayed to initialize because _nSH is set in derived class
  for (unsigned int j = 0; j < _npoints; ++j)
    _sample_point_data.push_back(MultiIndex<Real>({_I, _G, _nSH}));
}

void
NeutronicsSpectrumSamplerBase::initialize()
{
  meshChanged();
}

void
NeutronicsSpectrumSamplerBase::finalize()
{
  for (unsigned j = 0; j < _npoints; ++j)
    for (unsigned int i = 0; i < _I; ++i)
      for (unsigned int g = 0; g < _G; ++g)
        for (unsigned int p = 0; p < _nSH; ++p)
          gatherSum(_sample_point_data[j]({i, g, p}));

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
        for (unsigned int p = 0; p < _nSH; ++p)
          _sample_point_data[j]({i, g, p}) += uo._sample_point_data[j]({i, g, p});
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
  mooseError("Isotope name " << s << " cannot be converted with the localStringToZaid conversion method.");
  return 0;
}
