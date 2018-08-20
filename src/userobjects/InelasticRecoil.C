/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "InelasticRecoil.h"
#include "MagpieUtils.h"
#include "Function.h"

// gsl includes
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

registerMooseObject("MagpieApp", InelasticRecoil);

template <>
InputParameters
validParams<InelasticRecoil>()
{
  InputParameters params = validParams<ScatteringRecoilCrossSection>();
  params.addClassDescription(
      "Computes recoil cross sections for inelastic scattering events. Allows output to csv.");
  params.addRequiredParam<std::vector<Real>>("Q",
                                             "The Q values for all inelastic reaction channels.");
  return params;
}

InelasticRecoil::InelasticRecoil(const InputParameters & parameters)
  : ScatteringRecoilCrossSection(parameters),
    _q_values(getParam<std::vector<Real>>("Q")),
    _n_levels(_q_values.size())
{
  if (_n_levels != _scattering_cross_section.size())
    mooseError("Number of scattering_xs entries must be equal to number of Q values. Each "
               "inelastic mode is represented by one Q value.");
}

Real
InelasticRecoil::getLabCosine(Real E, Real T, Real Q) const
{
  Real mu_C = getCMCosine(E, T, Q);
  /*
   *  need to avoid NaN here from E < Eth , this can occur on input when
   * El (lower group bound) is supplied to compute LAB mu_min/max, note that
   * this routine returns one; this is consistent with the case Q = 0, E = 0 => mu_C = -1
   *
   * We also need to guard against the case mu_C > 1 which can occur when mu_L_max is computed
   * Both P1 = (E_l, T_l) and P2 = (E_l, T_u) can correspond to impermissible states; compare Fig. 3
   * in writeup. P1 and P2 can be outside of the permissible range of energies given by the green
   * and blue curves.
   */
  Real Eth = std::abs(Q) * (_atomic_mass + 1.0) / _atomic_mass;
  if (E <= Eth || mu_C >= 1)
    return 1;
  Real d = Eth / E;
  Real f = std::sqrt(1.0 - mu_C * mu_C) / (1.0 / std::sqrt(1 - d) - mu_C);
  return std::sqrt(1.0 / (1.0 + f * f));
}

Real
InelasticRecoil::getCMCosine(Real E, Real T, Real Q) const
{
  mooseAssert(E >= 0, "Negative neutron energy");
  mooseAssert(T >= 0, "Negative recoil energy");
  // need to avoid NaN here from E < Eth , this can occur on input when
  // El (lower group bound) is supplied to compute CM mu_min/max
  Real Eth = std::abs(Q) * (_atomic_mass + 1.0) / _atomic_mass;
  if (E <= Eth)
    return -1;
  Real d = Eth / E;
  Real mu = (1.0 - 2.0 * T / E / _gamma - 0.5 * d) / std::sqrt(1.0 - d);
  if (mu < -1)
    return -1;
  return mu;
}

Real
InelasticRecoil::getMaxRecoilEnergy(Real E, Real Et)
{
  if (E <= Et)
    return 0;
  Real delta = Et / E;
  return 0.5 * _gamma * E * (1.0 - 0.5 * delta + std::sqrt(1 - delta));
}

Real
InelasticRecoil::getMinRecoilEnergy(Real E, Real Et)
{
  if (E <= Et)
    return 0;
  Real delta = Et / E;
  return 0.5 * _gamma * E * (1.0 - 0.5 * delta - std::sqrt(1 - delta));
}

void
InelasticRecoil::execute()
{
  // inelastic scattering is isotropic in CM. To remind us that the scattering scattering law
  // is present, this variable is created and used
  const Real scattering_law = 0.5;

  for (unsigned int level = 0; level < _n_levels; ++level)
  {
    Real Q = _q_values[level];
    Real threshold_energy = std::abs(Q) * (_atomic_mass + 1.0) / _atomic_mass;

    // Loop over all Legendre orders to calculate the coefficients
    for (unsigned int l = 0; l < _L + 1; ++l)
    {
      // Loop over all the neutron energy groups
      for (unsigned int g = 0; g < _G; ++g)
      {
        // Gets the upper and lower energy values of that group
        Real E_u = _neutron_energy_limits[g];
        Real E_l = _neutron_energy_limits[g + 1];

        // we can skip this energy group of maximum energy is below or equal to threshold
        if (E_u <= threshold_energy)
          continue;

        // Loop over all the neutron energies within group g
        for (unsigned int i_E = 0; i_E < _quad_points.size(); ++i_E)
        {
          // Get the energies in Gaussian quadrature points and their respective weights
          Real E = 0.5 * (E_u - E_l) * _quad_points[i_E] + 0.5 * (E_u + E_l);
          Real w_E = 0.5 * _quad_weights[i_E] * (E_u - E_l);

          // This is the dimensionless neutron energy [comp. writeup]
          Real delta = threshold_energy / E;

          // Calculate maximum amount of energy transferred to the recoil atom
          Real T_min = getMinRecoilEnergy(E, threshold_energy);
          Real T_max = getMaxRecoilEnergy(E, threshold_energy);

          // Loop over all the possible recoil energy bins
          for (unsigned int t = 0; t < _T; ++t)
          {
            // Gets the upper and lower energy values of that bin
            Real T_u = _recoil_energy_limits[t];
            Real T_l = _recoil_energy_limits[t + 1];

            // Adjust the upper and lower limit according to inelastic scattering limits
            // [this is described in the document scattering_recoils.pdf]
            if (E_l < threshold_energy)
              E_l = threshold_energy;

            // both called with E_u is intendend [see scattering_recoils.pdf]
            T_u = std::min(T_u, getMaxRecoilEnergy(E_u, threshold_energy));
            T_l = std::max(T_l, getMinRecoilEnergy(E_u, threshold_energy));

            // if T_u <= T_l this recoil range is impossible
            if (T_u <= T_l)
              continue;

            // Compute the range of lab cosines. NOTE: this is more complicated than
            // for elastic scattering because mu_L assumes a minimum for T = |Q| / (A+1)
            Real mu_L_max = std::max(getLabCosine(E_l, T_l, Q), getLabCosine(E_l, T_u, Q));
            Real mu_L_min;
            Real Tp = std::abs(Q) / (_atomic_mass + 1);
            if (T_l <= Tp && T_u >= Tp)
              mu_L_min = std::sqrt(threshold_energy / E_u);
            else
            {
              if (T_u < Tp)
                mu_L_min = getLabCosine(E_u, T_u, Q);
              else if (T_l > Tp)
                mu_L_min = getLabCosine(E_u, T_l, Q);
              else
                mooseError("Should never get here");
            }

            // Save maximum and mininum lab frame cosine values
            _mu_L[t][g][0] = mu_L_min;
            _mu_L[t][g][1] = mu_L_max;

            for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
            {
              // Get the energies in Gaussian quadrature points and their respective weights
              Real T = 0.5 * (T_u - T_l) * _quad_points[i_T] + 0.5 * (T_u + T_l);

              // T must still adhere to the limits T_min and T_max!
              if (T <= T_min || T >= T_max)
                continue;

              Real w_T = 0.5 * _quad_weights[i_T] * (T_u - T_l);

              // Calculate cosine of recoil angle in the Lab frame according to geometry rules
              Real mu_L = getLabCosine(E, T, Q);

              /*
               * Calculate contribution to cross section coefficients
               * mu_L is scaled from its possible range of values [mu_L_min, mu_L_max] to fit the
               * interval [-1,1] of the Legendre polynomials
               */
              Real scaled_mu_L = 2 * (mu_L - mu_L_min) / (mu_L_max - mu_L_min) - 1;

              Real jacobian_determinant = 2.0 / (_gamma * E * std::sqrt(1 - delta));
              _recoil_coefficient[l][t][g] +=
                  1 / _xi_g[g] * _scattering_cross_section[level]->value(E, Point()) *
                  _neutron_spectrum.value(E, Point()) * scattering_law * jacobian_determinant *
                  gsl_sf_legendre_Pl(l, scaled_mu_L) * w_T * w_E;
            }
          }
        }
      }
    }
  }
}

#endif
