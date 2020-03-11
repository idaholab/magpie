/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef GSL_ENABLED

#include "DPAUserObjectBase.h"
#include "PolyatomicDisplacementFunction.h"

InputParameters
DPAUserObjectBase::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addParam<Real>("irradiation_time", "Irradiation time used ");
  MultiMooseEnum damage_reaction_types("elastic inelastic");
  params.addParam<MultiMooseEnum>("damage_reaction_types",
                                  damage_reaction_types,
                                  "The neutron reaction causing radiation damage");
  params.addParam<std::vector<Real>>("Z", "The atomic numbers of all present isotopes");
  params.addParam<std::vector<Real>>("A", "The mass numbers of all present isotopes");
  params.addParam<std::vector<Real>>("number_densities",
                                     "The number densities of all present isotopes");
  params.addParam<std::vector<Real>>("energy_group_boundaries",
                                     "The neutron flux energy group boundaries in units of eV "
                                     "starting with the highest energy group");
  params.addParam<std::vector<Real>>("scalar_flux",
                                     "The values of the neutron scalar flux by energy group "
                                     "starting with the highest energy group");
  params.addParam<std::vector<std::vector<Real>>>(
      "cross_section",
      "The values of the cross sections. Each semicolon separated vector contains cross sections "
      "for a particular nuclide and reaction type. Each entry must be number of energy groups "
      "entries long. One vector must be provided for each "
      "nuclide/reaction type combination. The ordering is as follows: if there are reaction types"
      "a and b, and nuclides i, k, and l, the ordering will be xs_ai; xs_ak; xs_al; xs_bi; xs_bk, "
      "xs_bl");
  params.addParam<std::vector<std::vector<Real>>>(
      "Q", "The Q values by reaction type and then isotope. Assumed zero if not provided.");
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

DPAUserObjectBase::DPAUserObjectBase(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _is_transient_irradiation(!isParamValid("irradiation_time")),
    _irradiation_time(_is_transient_irradiation ? 0 : getParam<Real>("irradiation_time")),
    _neutron_reaction_types(getParam<MultiMooseEnum>("damage_reaction_types")),
    _nr(_neutron_reaction_types.size())
{
  if (_nr == 0)
    paramError("damage_reaction_types", "At least one damage mechanism must be provided");

  // get parameters if they are provided
  if (isParamValid("Z"))
    _atomic_numbers = getParam<std::vector<Real>>("Z");
  if (isParamValid("A"))
    _mass_numbers = getParam<std::vector<Real>>("A");
  if (isParamValid("number_densities"))
    _number_densities = getParam<std::vector<Real>>("number_densities");
  if (isParamValid("energy_group_boundaries"))
    _energy_group_boundaries = getParam<std::vector<Real>>("energy_group_boundaries");
  if (isParamValid("scalar_flux"))
    _scalar_flux = getParam<std::vector<Real>>("scalar_flux");

  if (isParamValid("cross_section"))
  {
    std::vector<std::vector<Real>> xs = getParam<std::vector<std::vector<Real>>>("cross_section");

    // we know _nr so the length of cross sections should be divisible by _nr
    unsigned int n = xs.size();
    if (n % _nr != 0)
      paramError("cross_section",
                 "The number of entries in the cross_section param must be divisible by number of "
                 "reaction types ",
                 _nr);

    unsigned int n_nuc = n / _nr;
    _cross_sections.resize(_nr);
    unsigned int p = 0;
    for (unsigned int r = 0; r < _nr; ++r)
    {
      _cross_sections[r].resize(n_nuc);
      for (unsigned int j = 0; j < n_nuc; ++j)
      {
        _cross_sections[r][j] = xs[p];
        ++p;
      }
    }
  }
}

void
DPAUserObjectBase::prepare()
{
  // This object support dynamic changes of number densities
  // appearance and disappearance of isotopes, and change in
  // cross sections. Therefore, this function needs to check
  // and redo a lot of things that would usually be done in
  // the constructor.

  // checks on scalar flux & energy groups
  _ng = _scalar_flux.size();
  if (_ng == 0)
    mooseError("The number of provided scalar fluxes is zero. Provide using parameter scalar_flux "
               "or set programmatically");

  if (_energy_group_boundaries.size() != _ng + 1)
    mooseError("The size of the energy group boundary is ",
               _energy_group_boundaries.size(),
               " but it must have exactly number of energy groups plus one entry ",
               _ng + 1);

  for (unsigned int j = 1; j < _energy_group_boundaries.size(); ++j)
    if (_energy_group_boundaries[j - 1] <= _energy_group_boundaries[j])
      mooseError("Entries of energy_group_boundaries must be strictly decreasing");

  // check on A, Z, number densities
  initAZNHelper();

  // checks on cross sections
  if (_cross_sections.size() != _nr)
    mooseError("Leading dimension of _cross_sections is ",
               _cross_sections.size(),
               " but should be equal to number of damage reaction types ",
               _nr);
  for (unsigned int i = 0; i < _nr; ++i)
  {
    if (_cross_sections[i].size() != _atomic_numbers.size())
      mooseError("Second dimension of _cross_sections for index ",
                 i,
                 " does not match number of isotope ",
                 _atomic_numbers.size());
    for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
      if (_cross_sections[i][j].size() != _ng)
        mooseError("Third dimension of _cross_sections for indices ",
                   i,
                   " ",
                   j,
                   " has ",
                   _cross_sections[i][j].size(),
                   " entries when number of groups is ",
                   _ng);
  }

  // get Q values, need to do this here to be sure that number of isotopes is set
  if (isParamValid("Q"))
  {
    _q_values = getParam<std::vector<std::vector<Real>>>("Q");
    if (_q_values.size() != _nr)
      paramError("Q", "Leading dimension must be of the same size as reaction types.");
    for (unsigned int j = 0; j < _nr; ++j)
      if (_q_values[j].size() != _atomic_numbers.size())
        paramError("Q", "Second dimension must be of same size as isotopes.");
  }
  else
  {
    _q_values.resize(_nr);
    for (unsigned int j = 0; j < _nr; ++j)
      _q_values[j].resize(_atomic_numbers.size());
  }
}

void
DPAUserObjectBase::initAZNHelper()
{
  unsigned int ns = getAtomicNumbers().size();
  if (ns == 0)
    mooseError("The size of the atomic number array is zero. Set Z parameter or set atomic number "
               "programmatically");

  if (_mass_numbers.size() != ns)
    mooseError("Size of mass number array does not match size of atomic number array. Size of mass "
               "numbers ",
               _mass_numbers.size(),
               ", size of atomic numbers ",
               ns);

  if (_number_densities.size() != ns)
    mooseError("Size of mass number array does not match size of atomic number array. Size of mass "
               "numbers ",
               _number_densities.size(),
               ", size of atomic numbers ",
               ns);
}

std::vector<unsigned int>
DPAUserObjectBase::getAtomicNumbers() const
{
  // in this kind we need to work a bit since
  // VPPs are Reals and we need to return unsigned int
  std::vector<unsigned int> Zs(_atomic_numbers.size());
  for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
  {
    unsigned int i = static_cast<unsigned int>(_atomic_numbers[j]);
    if (std::abs(_atomic_numbers[j] - i) > 1e-12)
      mooseError("Entry ", j, ":", _atomic_numbers[j], " is not a non-negative integer.");
    Zs[j] = i;
  }
  return Zs;
}

std::vector<Real>
DPAUserObjectBase::getMassNumbers() const
{
  return _mass_numbers;
}

std::vector<Real>
DPAUserObjectBase::getNumberFractions() const
{
  // need to normalize here
  Real sum = 0;
  for (unsigned int j = 0; j < _number_densities.size(); ++j)
    sum += _number_densities[j];
  std::vector<Real> nf(_number_densities.size());
  for (unsigned int j = 0; j < _number_densities.size(); ++j)
    nf[j] = _number_densities[j] / sum;
  return nf;
}

Real
DPAUserObjectBase::getMaxEnergy() const
{
  // between elastic & inelastic scattering, elastic scattering
  // allows the largest recoil energy = gamma * E_neutron
  Real max_neutron_energy = _energy_group_boundaries[0];

  // find maximum gamma
  Real min_mass = std::numeric_limits<Real>::max();
  for (unsigned int j = 0; j < _mass_numbers.size(); ++j)
    if (min_mass > _mass_numbers[j])
      min_mass = _mass_numbers[j];
  Real gamma = 4 * min_mass / (min_mass + 1) / (min_mass + 1);
  return gamma * max_neutron_energy;
}

bool
DPAUserObjectBase::changed() const
{
  // the size could have changed
  if (_atomic_numbers.size() != _atomic_numbers_old.size() ||
      _mass_numbers.size() != _mass_numbers_old.size() ||
      _number_densities.size() != _number_densities_old.size())
    return true;

  // the size is still the same
  for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
    if (!MooseUtils::absoluteFuzzyEqual(_atomic_numbers[j], _atomic_numbers_old[j], _tol) ||
        !MooseUtils::absoluteFuzzyEqual(_mass_numbers[j], _mass_numbers_old[j], _tol) ||
        !MooseUtils::absoluteFuzzyEqual(_number_densities[j], _number_densities_old[j], _tol))
      return true;

  return false;
}

void
DPAUserObjectBase::accumulateDamage()
{
  if (changed())
  {
    initAZNHelper();
    onCompositionChanged();

    // copy over A, Z, N to A, Z, N_old
    _atomic_numbers_old.resize(_atomic_numbers.size());
    _mass_numbers_old.resize(_atomic_numbers.size());
    _number_densities_old.resize(_atomic_numbers.size());
    for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
    {
      _atomic_numbers_old[j] = _atomic_numbers[j];
      _mass_numbers_old[j] = _mass_numbers[j];
      _number_densities_old[j] = _number_densities[j];
    }
  }

  // compute total number density
  Real total_number_density = 0;
  for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
    total_number_density += _number_densities[j];

  // damage function have been computed
  // so now we can start computing dpa
  if (!_is_transient_irradiation)
  {
    _dpa = 0;
    _partial_dpa.clear();
  }

  Real del_t = _is_transient_irradiation ? _dt : _irradiation_time;

  for (unsigned int x = 0; x < _nr; ++x)
  {
    for (unsigned int g = 0; g < _ng; ++g)
    {
      Real del_E = _energy_group_boundaries[g] - _energy_group_boundaries[g + 1];

      for (unsigned int i = 0; i < _atomic_numbers.size(); ++i)
        for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
        {
          Real v = del_t / del_E * _scalar_flux[g] * _number_densities[i] / total_number_density *
                   _cross_sections[x][i][g] * neutronDamageEfficiency(i, j, g, x);

          // total dpa is incremented regardless
          _dpa += v;

          // partial dpa: if the entry exists, add to it, otherwise create it
          const auto & p = _partial_dpa.find(zaidHelper(_atomic_numbers[i], _mass_numbers[i]));
          if (p == _partial_dpa.end())
            _partial_dpa[zaidHelper(_atomic_numbers[i], _mass_numbers[i])] = v;
          else
            p->second += v;
        }
    }
  }
}

Real
DPAUserObjectBase::neutronDamageEfficiency(unsigned int i,
                                           unsigned int j,
                                           unsigned int g,
                                           unsigned int x) const
{
  std::vector<Real> points;
  std::vector<Real> weights;
  PolyatomicDisplacementFunction::gslQuadRule(100, points, weights);

  if (_neutron_reaction_types[x] == "elastic")
  {
    Real gamma = 4 * _mass_numbers[i] / (_mass_numbers[i] + 1) / (_mass_numbers[i] + 1);
    Real from_E = _energy_group_boundaries[g + 1];
    Real to_E = _energy_group_boundaries[g];
    Real delta_E = to_E - from_E;

    // integral over E
    Real integral = 0;
    for (unsigned int qp = 0; qp < points.size(); ++qp)
    {
      Real energy = 0.5 * (points[qp] + 1) * delta_E + from_E;
      Real w = 0.5 * weights[qp] * delta_E;
      integral += w * integralDamageFunction(gamma * energy, i, j) / energy;
    }

    return integral / gamma;
  }
  else if (_neutron_reaction_types[x] == "inelastic")
  {
    // warn if q value does not make sense
    if (_q_values[x][i] >= 0)
      mooseDoOnce(mooseWarning("Q value is greater or equal to zero for inelastic scattering."));

    Real d = std::abs(_q_values[x][i]) * (_mass_numbers[i] + 1) / _mass_numbers[i];
    Real gamma = 4 * _mass_numbers[i] / (_mass_numbers[i] + 1) / (_mass_numbers[i] + 1);
    Real from_E = _energy_group_boundaries[g + 1];
    Real to_E = _energy_group_boundaries[g];
    Real delta_E = to_E - from_E;

    // integral over E
    Real integral = 0;
    for (unsigned int qp = 0; qp < points.size(); ++qp)
    {
      Real energy = 0.5 * (points[qp] + 1) * delta_E + from_E;
      Real Delta = d / energy;

      // make sure threshold energy is honored
      if (Delta > 1)
        continue;

      Real sr = std::sqrt(1 - Delta);
      Real Tmin = gamma / 2 * energy * (1 - sr - 0.5 * Delta);
      Real Tmax = gamma / 2 * energy * (1 + sr - 0.5 * Delta);
      Real w = 0.5 * weights[qp] * delta_E;
      integral += w * (integralDamageFunction(Tmax, i, j) - integralDamageFunction(Tmin, i, j)) /
                  (energy * sr);
    }

    return integral / gamma;
  }

  mooseError("Neutron reaction type not recognized. Should never get here.");
  return 0;
}

std::string
DPAUserObjectBase::zaidHelper(unsigned int Z, Real A) const
{
  std::stringstream ss;
  ss << std::setprecision(4) << Z << A;
  return ss.str();
}

Real
DPAUserObjectBase::getPartialDPA(unsigned int Z, Real A) const
{
  const auto & p = _partial_dpa.find(zaidHelper(Z, A));
  if (p == _partial_dpa.end())
    return 0;
  return p->second;
}

#endif
