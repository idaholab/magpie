/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef GSL_ENABLED
#ifdef RATTLESNAKE_ENABLED

#include "NeutronicsSpectrumSamplerSN.h"
#include "IsotopeUtilities.h"
#include "ElasticRecoil.h"

// gsl includes
#include "gsl/gsl_sf_legendre.h"

#include <algorithm>

registerMooseObject("MagpieApp", NeutronicsSpectrumSamplerSN);

InputParameters
NeutronicsSpectrumSamplerSN::validParams()
{
  InputParameters params = NeutronicsSpectrumSamplerBase::validParams();
  params.addRequiredCoupledVar("angular_variables",
                               "Angular fluxes, dimension G x M (# angular directions).");
  params.addRequiredParam<std::vector<std::string>>("recoil_isotope_names",
                                                    "The list of recoil isotope names e.g. U235.");
  params.addRequiredParam<UserObjectName>("aqdata", "Angular quadrature user data.");
  params.addRequiredParam<std::vector<UserObjectName>>(
      "recoil_cross_sections", "Recoil cross section UserObject names. Size = npoints x nisotopes");
  params.addParam<unsigned int>(
      "nmu", 10, "Number of polar subdivisions for angular dependence of PDF");
  params.addParam<unsigned int>(
      "nphi", 10, "Number of azimuthal subdivisions for angular dependence of PDF");
  params.addClassDescription("Computes PDFs for reactions except fission that can be used for "
                             "sampling PKAs in coupled BCMC simulations.\n User match match target "
                             "isotopes with recoil isotopes and provide transfer-like recoil cross "
                             "section data.");
  return params;
}

NeutronicsSpectrumSamplerSN::NeutronicsSpectrumSamplerSN(const InputParameters & parameters)
  : NeutronicsSpectrumSamplerBase(parameters),
    _recoil_isotope_names(getParam<std::vector<std::string>>("recoil_isotope_names")),
    _aq(getUserObject<AQData>("aqdata").aq()),
    _ndir(_aq.getNQuadratures()),
    _shm(SHCoefficients(_aq, _L))
{
  _nmu = getParam<unsigned int>("nmu");
  _nphi = getParam<unsigned int>("nphi");

  // check recoil cross section length
  std::vector<UserObjectName> names =
      getParam<std::vector<UserObjectName>>("recoil_cross_sections");
  if (names.size() != _npoints * _I)
    mooseError(
        "recoil_cross_sections must be a vector of length npoints x nisotopes of UserObjectNames");

  // check recoil ZAIDs
  if (_recoil_isotope_names.size() != _I)
    mooseError("recoil_isotope_names have wrong length.");
  unsigned int Z, A;
  for (unsigned int i = 0; i < _I; ++i)
    YAKXS::Utility::getAZ(_recoil_isotope_names[i], A, Z);

  // get angular fluxes
  _angular_flux.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
  {
    _angular_flux[g].resize(_ndir);
    for (unsigned int dir = 0; dir < _ndir; ++dir)
      _angular_flux[g][dir] = &coupledValue("angular_variables", g * _ndir + dir);
  }

  // allocate and assign recoil cross section
  unsigned int p = 0;
  _recoil_cross_sections.resize(_npoints);
  for (unsigned j = 0; j < _npoints; ++j)
  {
    _recoil_cross_sections[j].resize(_I);
    for (unsigned int i = 0; i < _I; ++i)
    {
      _recoil_cross_sections[j][i] = &getUserObjectByName<ElasticRecoil>(names[p]);
      ++p;
    }
  }
}

Real
NeutronicsSpectrumSamplerSN::totalRecoilRate(unsigned int point_id,
                                             const std::string & target_isotope) const
{
  Real rate = 0.0;
  Real dmu = 2.0 / Real(_nmu);
  Real dphi = 2.0 * libMesh::pi / Real(_nphi);
  auto it = std::find(_target_isotope_names.begin(), _target_isotope_names.end(), target_isotope);
  if (it == _target_isotope_names.end())
    mooseError("Isotope ", target_isotope, "does not exist");
  unsigned int target_isotope_id = std::distance(_target_isotope_names.begin(), it);
  for (unsigned int g = 0; g < _G; ++g)
    for (unsigned int p = 0; p < _nmu; ++p)
      for (unsigned int q = 0; q < _nphi; ++q)
        rate += dmu * dphi * _sample_point_data[point_id]({target_isotope_id, g, p, q});
  return rate;
}

Real
NeutronicsSpectrumSamplerSN::computeRadiationDamagePDF(unsigned int i,
                                                       unsigned int g,
                                                       unsigned int p,
                                                       unsigned int q)
{
  Real a = 0.0;

  // find the mu and phi ranges in the current bin
  Real lower_phi = p * 2.0 * libMesh::pi / Real(_nphi);
  Real upper_phi = (p + 1) * 2.0 * libMesh::pi / Real(_nphi);
  Real lower_mu = q * 2.0 / Real(_nmu) - 1.0;
  Real upper_mu = (q + 1) * 2.0 / Real(_nmu) - 1.0;

  // TODO for better accuracy we may use a quadrature rule over the intervals;
  // for now the midpoint is good enough
  Real mu = 0.5 * (lower_mu + upper_mu);
  Real phi = 0.5 * (lower_phi + upper_phi);
  RealVectorValue omega_T(
      mu, std::cos(phi) * std::sqrt(1 - mu * mu), std::sin(phi) * std::sqrt(1 - mu * mu));

  for (unsigned int gp = 0; gp < _G; ++gp)
    for (unsigned int dir = 0; dir < _ndir; ++dir)
    {
      Real mu_lab = omega_T * _aq.getDirectionInRV(dir);
      Real mu_lab_min = _recoil_cross_sections[_current_point][i]->getMinRecoilCosine(gp, g);
      Real mu_lab_max = _recoil_cross_sections[_current_point][i]->getMaxRecoilCosine(gp, g);

      // check if mu_lab is permissible
      if (mu_lab > mu_lab_max || mu_lab < mu_lab_min)
        continue;

      Real mu_transformed = 2.0 * (mu_lab - mu_lab_min) / (mu_lab_max - mu_lab_min) - 1.0;
      unsigned int L_data = _recoil_cross_sections[_current_point][i]->legendreOrder();

      if (L_data < _L)
        mooseDoOnce(mooseWarning("L is larger than the legendre order of the provided data"));

      for (unsigned int l = 0; l <= std::min(_L, L_data); ++l)
        a += (*_number_densities[i])[_qp] * _aq.getWeight(dir) *
             gsl_sf_legendre_Pl(l, mu_transformed) *
             _recoil_cross_sections[_current_point][i]->getSigmaRecoil(gp, g, l) *
             (*_angular_flux[gp][dir])[_qp];
    }
  return a;
}

#endif // RATTLESNAKE_ENABLED
#endif // GSL_ENABLED
