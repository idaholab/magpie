/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "ElasticRecoil.h"
#include "MagpieUtils.h"
#include "Function.h"

// gsl includes
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

registerMooseObject("MagpieApp", ElasticRecoil);

template <>
InputParameters
validParams<ElasticRecoil>()
{
  InputParameters params = validParams<ScatteringRecoilCrossSection>();
  params.addClassDescription(
      "Computes recoil cross sections for elastic scattering events. Allows output to csv.");
  params.addRequiredParam<FunctionName>("scattering_law",
                                        "Function representing the scattering law for neutrons");
  return params;
}

ElasticRecoil::ElasticRecoil(const InputParameters & parameters)
  : ScatteringRecoilCrossSection(parameters), _scattering_law(getFunction("scattering_law"))
{
  if (_scattering_cross_section.size() != 1)
    mooseError("ElasticRecoil only allows a single input for scattering_xs parameter. Multiple "
               "entries are reserved for inelastic scattering modes.");
}

Real
ElasticRecoil::getLabCosine(Real E, Real T, Real /*Q*/) const
{
  Real mu_C = getCMCosine(E, T);
  return std::sqrt((1 - mu_C) / 2);
}

Real
ElasticRecoil::getCMCosine(Real E, Real T, Real /*Q*/) const
{
  mooseAssert(E >= 0, "Negative neutron energy");
  mooseAssert(T >= 0, "Negative recoil energy");
  // just need to avoid NaN here, E=0 can occur on input when lower energy
  // bin boundary is E = 0
  if (E == 0)
    return -1;
  // ensure proper boundaries for CM cosine
  Real mu = 1 - 2 * T / E / _gamma;
  if (mu < -1)
    return -1;
  return mu;
}

void
ElasticRecoil::execute()
{
  // Loop over all Legendre orders to calculate the coefficients
  for (unsigned int l = 0; l < _L + 1; ++l)
  {
    // Loop over all the neutron energy groups
    for (unsigned int g = 0; g < _G; ++g)
    {
      // Gets the upper and lower energy values of that group
      Real E_u = _neutron_energy_limits[g];
      Real E_l = _neutron_energy_limits[g + 1];

      // Loop over all the neutron energies within group g
      for (unsigned int i_E = 0; i_E < _quad_points.size(); ++i_E)
      {
        // Get the energies in Gaussian quadrature points and their respective weights
        Real E = 0.5 * (E_u - E_l) * _quad_points[i_E] + 0.5 * (E_u + E_l);
        Real w_E = 0.5 * _quad_weights[i_E] * (E_u - E_l);

        // Calculate maximum amount of energy transferred to the recoil atom
        Real T_max = _gamma * E;

        // Loop over all the possible recoil energy bins
        for (unsigned int t = 0; t < _T; ++t)
        {
          // Gets the upper and lower energy values of that bin
          Real T_u = _recoil_energy_limits[t];
          Real T_l = _recoil_energy_limits[t + 1];

          /*
           * Calculate possible range of angles according to neutron energy group
           * and recoil energy bin. This approach avoids unphysical behavior (negative
           * values due to convergence issue of expansion) of recoil cross section.
           * The subscript L corresponds to the Lab frame.
           */
          Real mu_L_min = getLabCosine(E_u, T_l);
          Real mu_L_max = getLabCosine(E_l, T_u);

          // Save maximum and mininum Lab frame cosine values
          _mu_L[t][g][0] = mu_L_min;
          _mu_L[t][g][1] = mu_L_max;

          /*
           * Elastic scaterring case III: T_max < T_l, stop
           * When the maximum recoil energy is lower than the lowest energy of the recoil energy bin
           * t
           */
          if (T_max < T_l)
            continue;

          /*
           * Elastic scaterring case II: T_max inside T bin, refit
           * When the maximum recoil energy is whitin the recoil energy bin t, we need to refit
           * the quadrature rule for the new interval between T_l and T_max
           */
          if (T_max > T_l && T_max < T_u)
            T_u = T_max;

          /*
           * Elastic scaterring case I: T_max > interval, ok
           * When the maximum recoil energy is greater than the highest energy of the recoil energy
           * bin t Case II also utilizes this piece of code with T_u = T_max
           */
          for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
          {
            // Get the energies in Gaussian quadrature points and their respective weights
            Real T = 0.5 * (T_u - T_l) * _quad_points[i_T] + 0.5 * (T_u + T_l);

            Real w_T = 0.5 * _quad_weights[i_T] * (T_u - T_l);

            // Calculate cosine of recoil angle in the CM frame given the neutron and recoil atom
            // energies
            Real mu_C = getCMCosine(E, T);

            // Calculate cosine of recoil angle in the Lab frame according to geometry rules
            Real mu_L = getLabCosine(E, T);

            /*
             * Calculate contribution to cross section coefficients
             * mu_L is scaled from its possible range of values [mu_L_min, mu_L_max] to fit the
             * interval [-1,1] of the Legendre polynomials
             */
            Real scaled_mu_L = 2 * (mu_L - mu_L_min) / (mu_L_max - mu_L_min) - 1;
            _recoil_coefficient[l][t][g] +=
                1 / _xi_g[g] * _scattering_cross_section[0]->value(E, Point()) *
                _neutron_spectrum.value(E, Point()) * _scattering_law.value(mu_C, Point()) * 2.0 /
                _gamma / E * gsl_sf_legendre_Pl(l, scaled_mu_L) * w_T * w_E;
          }
        }
      }
    }
  }
}

#endif
