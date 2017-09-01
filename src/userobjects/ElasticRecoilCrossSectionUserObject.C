
#include "ElasticRecoilCrossSectionUserObject.h"
#include "MagpieUtils.h"
#include "Function.h"

// libmesh includes
#include "libmesh/quadrature.h"

// gsl includes
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

// Output includes
#include <cassert>
#include <fstream>
#include <stdexcept>

template <>
InputParameters
validParams<ElasticRecoilCrossSectionUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription("Calculate the recoil atom elastic cross section given the isotope and energy groups."
      "It outputs the coefficients for the Legendre expansion of the cross section according to the order of Legendre"
      "polynomials l, neutron energy groups g and recoil energy bins t. It also outputs the maximum and mininum"
      "cosines of the recoil atom scatter angle in the Laboratory frame.");
  params.addParam<unsigned int>(
      "quadrature_order", 400, "Quadrature order for integration");
  params.addParam<Real>(
      "atomic_mass", 1, "Atomic Mass of the isotope. Default to Hydrogen A = 1");
  params.addRequiredParam<std::vector<Real>>(
      "neutron_energy_limits", "Energy limits of the incident neutron in [eV] and descending order");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_limits", "Energy limits of the recoil atom in [eV] and descending order");
  params.addRequiredParam<FunctionName>(
      "neutron_spectrum","Function representing the reactor neutron spectrum");
  params.addRequiredParam<FunctionName>(
      "scattering_law", "Function representing the scattering law for neutrons");
  params.addRequiredParam<FunctionName>(
      "elastic_xs", "Function representing the neutron elastic cross section");
  params.addParam<unsigned int>(
      "legendre_order", 10, "Order of Legendre polynomials where n = 0, ..., 10. Default to P10");
  params.addParam<std::string>(
      "erxs_output_file_name", "Name of the output file with the cross section coefficients (.csv)");
  params.addParam<std::string>(
      "mu_L_output_file_name", "Name of the output file with the mininum and maximum mu_L (.csv)");
  return params;
}

ElasticRecoilCrossSectionUserObject::ElasticRecoilCrossSectionUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _csv_tolerance(1.0e-10),
    _quad_order(getParam<unsigned int>("quadrature_order")),
    _neutron_spectrum(getFunction("neutron_spectrum")),
    _scattering_law(getFunction("scattering_law")),
    _elastic_xs(getFunction("elastic_xs")),
    _L(getParam<unsigned int>("legendre_order")),
    _atomic_mass(getParam<Real>("atomic_mass")),
    _neutron_energy_limits(getParam<std::vector<Real>>("neutron_energy_limits")),
    _recoil_energy_limits(getParam<std::vector<Real>>("recoil_energy_limits"))
{
  if (isParamValid("erxs_output_file_name") ^ isParamValid("mu_L_output_file_name"))
    mooseError("erxs_output_file_name and mu_L_output_file_name must either both be present or absent");

  const gsl_integration_fixed_type * qtpe = gsl_integration_fixed_legendre;
  gsl_integration_fixed_workspace * workspace = gsl_integration_fixed_alloc(qtpe, _quad_order, -1.0, 1.0, 0.0, 0.0);
  for (unsigned int j = 0; j < _quad_order; ++j)
  {
    _quad_points.push_back(workspace->x[j]);
    _quad_weights.push_back(workspace->weights[j]);
  }
  gsl_integration_fixed_free(workspace);

  _alpha = std::pow(((_atomic_mass - 1) / (_atomic_mass + 1)), 2);
  _gamma = 4 * _atomic_mass / std::pow(( _atomic_mass + 1), 2);
  _G = _neutron_energy_limits.size() - 1;
  _T = _recoil_energy_limits.size() - 1;
}

void
ElasticRecoilCrossSectionUserObject::initialize()
{
  // Calculate neutron spectrum over group g
  _xi_g.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
  {
    Real E_u = _neutron_energy_limits[g];
    Real E_l = _neutron_energy_limits[g + 1];

    for (unsigned int p = 0; p < _quad_points.size(); ++p)
    {
      Real E = 0.5 * (E_u - E_l) * _quad_points[p] + 0.5 * (E_u + E_l);
      Real w_E = 0.5 * _quad_weights[p] * (E_u - E_l);
      _xi_g[g] += w_E * _neutron_spectrum.value(E, Point());
    }
  }
}

// Find the neutron energy group given a neutron energy
unsigned int
ElasticRecoilCrossSectionUserObject::findNeutronEnergyGroup(Real energy)
{
  for (unsigned int g = 0; g < _G; ++g)
  {
    if (energy < _neutron_energy_limits[g] && energy > _neutron_energy_limits[g + 1])
      return g;
  }
  mooseError("Should never get here");
}

void
ElasticRecoilCrossSectionUserObject::execute()
{
  // Size the cross section array
  _erxs_coeff.resize(_L + 1);
  for (unsigned int l = 0; l < _L + 1; ++l)
  {
    _erxs_coeff[l].resize(_T);
    for (unsigned int t = 0; t < _T; ++t)
      _erxs_coeff[l][t].resize(_G);
  }

  // Size the lab frame cosine array
  _save_mu_L.resize(_T);
  for (unsigned int t = 0; t < _T; ++t)
  {
    _save_mu_L[t].resize(_G);
    for (unsigned int g = 0; g < _G; ++g)
      _save_mu_L[t][g].resize(2);
  }

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
           * The subscript C means the Center of Mass (CM) frame and L, the Lab frame.
           */
          Real mu_C_max = 1 - 2 * T_u / E_l / _gamma;
          Real mu_C_min = 1 - 2 * T_l / E_u / _gamma;
          // both max and min can be < -1 leading to a mu_lab > 1; this looks bad
          // in CSV file so fix it here
          if (mu_C_max < -1)
            mu_C_max = -1;
          if (mu_C_min < -1)
            mu_C_min = -1;
          Real mu_L_max = std::sqrt((1 - mu_C_max) / 2);
          Real mu_L_min = std::sqrt((1 - mu_C_min) / 2);

          // Save maximum and mininum Lab frame cosine values
          _save_mu_L[t][g][0] = mu_L_min;
          _save_mu_L[t][g][1] = mu_L_max;

          /*
           * Elastic scaterring case III: T_max < T_l, stop
           * When the maximum recoil energy is lower than the lowest energy of the recoil energy bin t
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
           * When the maximum recoil energy is greater than the highest energy of the recoil energy bin t
           * Case II also utilizes this piece of code with T_u = T_max
           */
          for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
          {
            // Get the energies in Gaussian quadrature points and their respective weights
            Real T = 0.5 * (T_u - T_l) * _quad_points[i_T] + 0.5 * (T_u + T_l);
            Real w_T = 0.5 * _quad_weights[i_T] * (T_u - T_l);

            // Calculate cosine of recoil angle in the CM frame given the neutron and recoil atom energies
            Real mu_C = 1 - 2 * T / (E * _gamma);

            // Calculate cosine of recoil angle in the Lab frame according to geometry rules
            Real mu_L = sqrt((1 - mu_C) / 2);

            /*
             * Calculate contribution to cross section coefficients
             * mu_L is scaled from its possible range of values [mu_L_min, mu_L_max] to fit the interval [-1,1]
             * of the Legendre polynomials
             */
             Real scaled_mu_L = 2 * (mu_L - mu_L_min) / (mu_L_max - mu_L_min) - 1;
            _erxs_coeff[l][t][g] += 1 / _xi_g[g] *
                                    _elastic_xs.value(E,Point()) *
                                    _neutron_spectrum.value(E, Point()) *
                                    _scattering_law.value(mu_C, Point()) *
                                    2.0 / _gamma / E *
                                    gsl_sf_legendre_Pl(l, scaled_mu_L) * w_T * w_E;
          }
        }
      }
    }
  }
}

void
ElasticRecoilCrossSectionUserObject::finalize()
{
  if (!(isParamValid("erxs_output_file_name") && isParamValid("mu_L_output_file_name")))
    return;
  _erxs_file_name = getParam<std::string>("erxs_output_file_name");
  _mu_L_file_name = getParam<std::string>("mu_L_output_file_name");

  /*
   * Write the elastic recoil cross section (g -> t) output file
   *
   *            t = 0              t = 1              t = ...
   *       l = 0, 1, 2, ...   l = 0, 1, 2, ...   l = 0, 1, 2, ...
   *  g 0
   *    1
   *    2
   *    .
   *    .
   *    .
   */
  std::ofstream output_file;
  output_file.open(_erxs_file_name);

  // print header
  output_file << "Neutron group,";
  for (unsigned int t = 0; t < _T; ++t)
    for (unsigned int l = 0; l < _L + 1; ++l)
    {
      output_file << "Recoil group t = " << t << " Moment order l = " << l;
      if (!(t == _T - 1 && l == _L))
        output_file << ",";
    }
  output_file << std::endl;

  // print data
  for (unsigned int g = 0; g < _G; ++g)
  {
    output_file << g << ",";
    for (unsigned int t = 0; t < _T; ++t)
      for (unsigned int l = 0; l < _L + 1; ++l)
      {
        output_file << csvPrint(_erxs_coeff[l][t][g]);
        if (!(t == _T - 1 && l == _L))
          output_file << ',';
      }
    output_file << std::endl;
  }
  output_file.close();

  /*
   * Writes output file with maximum and mininum cosines in the Lab frame (mu_L)
   * It follows same structure as ERXS outfile file, but saves the mu_L_min and
   * mu_L_max per g -> t combination.
   */
  std::ofstream output_file2;
  output_file2.open(_mu_L_file_name);
  // print header
  output_file2 << "Neutron group,";
  for (unsigned int t = 0; t < _T; ++t)
  {
    output_file2 << "Recoil group t = " << t << " mu_min," << "Recoil group t = " << t << " mu_max";
    if (t != _T - 1)
      output_file2 << ",";
  }
  output_file2 << std::endl;

  // print data
  for (unsigned int g = 0; g < _G; ++g)
  {
    output_file2 << g << ",";
    for (unsigned int t = 0; t < _T; ++t)
    {
      output_file2 << csvPrint(_save_mu_L[t][g][0]) << ',' << csvPrint(_save_mu_L[t][g][1]);
      if (t != _T - 1)
        output_file2 << ',';
    }
    output_file2 << std::endl;
  }
  output_file2.close();
}

Real
ElasticRecoilCrossSectionUserObject::csvPrint(Real value) const
{
  if (std::abs(value) < _csv_tolerance)
    return 0;
  return value;
}

Real
ElasticRecoilCrossSectionUserObject::getSigmaRecoil(unsigned int g, unsigned int t, unsigned int l) const
{
  mooseAssert(g < _G, "g is larger than _G [indexed from g = 0, ..., G - 1]");
  mooseAssert(t < _T, "t is larger than _T [indexed from t = 0, ..., T - 1]");
  if (l <= _L)
    return _erxs_coeff[l][t][g];
  return 0.0;
}

Real
ElasticRecoilCrossSectionUserObject::getMaxRecoilCosine(unsigned int g, unsigned t) const
{
  mooseAssert(g < _G, "g is larger than _G [indexed from g = 0, ..., G - 1]");
  mooseAssert(t < _T, "t is larger than _T [indexed from t = 0, ..., T - 1]");
  return _save_mu_L[t][g][0];
}

Real
ElasticRecoilCrossSectionUserObject::getMinRecoilCosine(unsigned int g, unsigned t) const
{
  mooseAssert(g < _G, "g is larger than _G [indexed from g = 0, ..., G - 1]");
  mooseAssert(t < _T, "t is larger than _T [indexed from t = 0, ..., T - 1]");
  return _save_mu_L[t][g][1];
}
