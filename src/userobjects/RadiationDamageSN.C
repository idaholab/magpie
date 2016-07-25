/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifdef RATTLESNAKE_ENABLED
#include "RadiationDamageSN.h"
#include "YakxsUtilities.h"

template<>
InputParameters validParams<RadiationDamageSN>()
{
  InputParameters params = validParams<RadiationDamageBase>();
  params.addRequiredCoupledVar("angular_variables", "Angular fluxes, dimension G x M (# angular directions).");
  params.addRequiredParam<std::vector<std::string> >("recoil_isotope_names", "The list of recoil isotope names e.g. U235.");
  // FIXME: This is not the permanent solution for providing recoil cross sections. It must be implemented as material
  // property.
  params.addRequiredParam<UserObjectName>("aqdata", "Angular quadrature user data.");
  params.addRequiredParam<std::vector<Real> >("recoil_cross_sections", "Recoil cross sections. Size = npoints x nisotopes x G x G x (L+1).");
  params.addClassDescription("Computes PDFs for reactions except fission that can be used for sampling PKAs in coupled BCMC simulations.\n User match match target isotopes with recoil isotopes and provide transfer-like recoil cross section data.");
  return params;
}

RadiationDamageSN::RadiationDamageSN(const InputParameters & parameters) :
    RadiationDamageBase(parameters),
    _recoil_isotope_names(getParam<std::vector<std::string> >("recoil_isotope_names")),
    _aq(getUserObject<AQData>("aqdata").aq()),
    _ndir(_aq.getNQuadratures()),
    _shm(SHCoefficients(_aq, _L)),
    _nSH = _shm.getNSH(); // Set the number of spherical harmonics
{
  // check recoil cross section length
  std::vector<Real> rxs = getParam<std::vector<Real> >("recoil_cross_sections");
  if (rxs.size() != _npoints * _I * _G * _G * (_L + 1))
    mooseError("recoil cross sections must be of length npoints x nisotopes x G**2 x (L+1)");

  // check recoil ZAIDs
  if (_recoil_isotope_names.size() != _I)
    mooseError("recoil_isotope_names have wrong length.");
  unsigned int Z, A;
  for (unsigned int i = 0; i < _I; ++i)
    YAKXS::Utility::getAZFromIsotopeName(_recoil_isotope_names[i], A, Z);

  // get angular fluxes
  _angular_flux.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
  {
    _angular_flux[g].resize(_ndir);
    for (unsigned int dir = 0; dir < _ndir; ++dir)
      _angular_flux[g][dir] = & coupledValue("angular_variables", g * _ndir + dir);
  }

  // allocate flux moments
  _flux_moment.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
    _flux_moment[g].resize(_nSH);

  // allocate and assign recoil cross section
  unsigned int p = 0;
  _recoil_cross_section.resize(_npoints);
  for (unsigned j = 0; j < _npoints; ++j)
  {
    _recoil_cross_section[j].resize(_I);
    for (unsigned int i = 0; i < _I; ++i)
    {
      _recoil_cross_section[j][i].resize(_L + 1);
      for (unsigned l = 0; l < _L + 1; ++l)
      {
        // first allocate the whole block
        _recoil_cross_section[j][i][l].resize(_G);
        for (unsigned int g = 0; g < _G; ++g)
          _recoil_cross_section[j][i][l][g].resize(_G);
        // set the values
        for (unsigned int g = 0; g < _G; ++g)
          for (unsigned int gp = 0; gp < _G; ++gp)
            _recoil_cross_section[j][i][l][gp][g] = rxs[p++];
      }
    }
  }
}

void
RadiationDamageSN::preComputeRadiationDamagePDF()
{
  // compute _flux_moments for current _qp
  for (unsigned int g = 0; g < _G; ++g)
    for (unsigned int p = 0; p < _nSH; ++p)
    {
       _flux_moment[g][p] = 0.0;
       for (unsigned int dir = 0; dir < _ndir; ++dir)
         _flux_moment[g][p] += _aq.getWeight(dir) * _shm.getSH(p, dir) * (*_angular_flux[g][dir])[_qp];
    }
}

Real
RadiationDamageSN::computeRadiationDamagePDF(unsigned int i, unsigned int g, unsigned int p)
{
  Real a = 0.0;
  for (unsigned int gp = 0; gp < _G; ++gp)
    a += _flux_moment[gp][p] * (*_number_densities[i])[_qp] * _recoil_cross_section[_current_point][i][_shm.p2l(p)][gp][g];
  return a;
}

#endif //RATTLESNAKE_ENABLED
