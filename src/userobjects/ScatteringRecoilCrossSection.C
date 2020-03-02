/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "ScatteringRecoilCrossSection.h"
#include "MagpieUtils.h"
#include "Function.h"

// gsl includes
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"

InputParameters
ScatteringRecoilCrossSection::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Base class for computing recoil scattering cross sections for a single isotope."
      "It outputs the coefficients for the Legendre expansion of the cross section up to the "
      "specified Legendre expansion order"
      "polynomials l, neutron energy groups g and recoil energy bins t. It also outputs the "
      "maximum and mininum"
      "cosines of the recoil atom in the laboratory frame.");
  params.addParam<unsigned int>("quadrature_order", 400, "Quadrature order for integration");
  params.addParam<Real>("atomic_mass", 1, "Atomic Mass of the isotope. Default to Hydrogen A = 1");
  params.addRequiredParam<std::vector<Real>>(
      "neutron_energy_limits",
      "Energy limits of the incident neutron in [eV] and descending order");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_limits", "Energy limits of the recoil atom in [eV] and descending order");
  params.addRequiredParam<FunctionName>("neutron_spectrum",
                                        "Function representing the reactor neutron spectrum");
  params.addRequiredParam<std::vector<FunctionName>>("scattering_xs",
                                                     "Functions representing the neutron cross "
                                                     "sections. NOTE: this is a vector to allow "
                                                     "separate inputs for inelastic modes.");
  params.addParam<unsigned int>(
      "legendre_order", 10, "Order of Legendre polynomials where n = 0, ..., 10. Default to P10");
  params.addParam<std::string>(
      "cross_section_output_filename",
      "Name of the output file with the cross section coefficients (.csv)");
  params.addParam<std::string>("mu_L_output_filename",
                               "Name of the output file with the mininum and maximum mu_L (.csv)");
  return params;
}

ScatteringRecoilCrossSection::ScatteringRecoilCrossSection(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _csv_tolerance(1.0e-10),
    _quad_order(getParam<unsigned int>("quadrature_order")),
    _neutron_spectrum(getFunction("neutron_spectrum")),
    _L(getParam<unsigned int>("legendre_order")),
    _atomic_mass(getParam<Real>("atomic_mass")),
    _neutron_energy_limits(getParam<std::vector<Real>>("neutron_energy_limits")),
    _recoil_energy_limits(getParam<std::vector<Real>>("recoil_energy_limits"))
{
  if (isParamValid("cross_section_output_filename") ^ isParamValid("mu_L_output_filename"))
    mooseError("cross_section_output_filename and mu_L_output_filename must either both be present "
               "or absent");

  // get scattering cross sections
  std::vector<FunctionName> xs_names = getParam<std::vector<FunctionName>>("scattering_xs");
  _scattering_cross_section.resize(xs_names.size());
  for (unsigned int j = 0; j < xs_names.size(); ++j)
    _scattering_cross_section[j] = &getFunctionByName(xs_names[j]);

  // set up integration rule
  auto * qp_table = gsl_integration_glfixed_table_alloc(_quad_order);
  _quad_points.resize(_quad_order);
  _quad_weights.resize(_quad_order);
  for (std::size_t j = 0; j < _quad_order; ++j)
  {
    double point, weight;
    gsl_integration_glfixed_point(-1.0, 1.0, j, &point, &weight, qp_table);
    _quad_points[j] = point;
    _quad_weights[j] = weight;
  }
  gsl_integration_glfixed_table_free(qp_table);

  _alpha = std::pow(((_atomic_mass - 1) / (_atomic_mass + 1)), 2);
  _gamma = 4 * _atomic_mass / std::pow((_atomic_mass + 1), 2);
  _G = _neutron_energy_limits.size() - 1;
  _T = _recoil_energy_limits.size() - 1;
}

void
ScatteringRecoilCrossSection::initialize()
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

  // Size the cross section array
  _recoil_coefficient.resize(_L + 1);
  for (unsigned int l = 0; l < _L + 1; ++l)
  {
    _recoil_coefficient[l].resize(_T);
    for (unsigned int t = 0; t < _T; ++t)
      _recoil_coefficient[l][t].assign(_G, 0.0);
  }

  // Size the lab frame cosine array and initialize to -1, both lab cosines -1 means
  // this g->t combination is impossible
  _mu_L.resize(_T);
  for (unsigned int t = 0; t < _T; ++t)
  {
    _mu_L[t].resize(_G);
    for (unsigned int g = 0; g < _G; ++g)
      _mu_L[t][g].assign(2, -1.0);
  }
}

void
ScatteringRecoilCrossSection::finalize()
{
  if (!(isParamValid("cross_section_output_filename") && isParamValid("mu_L_output_filename")))
    return;
  _recoil_xs_file_name = getParam<std::string>("cross_section_output_filename");
  _mu_L_file_name = getParam<std::string>("mu_L_output_filename");

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
  output_file.open(_recoil_xs_file_name);

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
        output_file << csvPrint(_recoil_coefficient[l][t][g]);
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
    output_file2 << "Recoil group t = " << t << " mu_min,"
                 << "Recoil group t = " << t << " mu_max";
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
      output_file2 << csvPrint(_mu_L[t][g][0]) << ',' << csvPrint(_mu_L[t][g][1]);
      if (t != _T - 1)
        output_file2 << ',';
    }
    output_file2 << std::endl;
  }
  output_file2.close();
}

// Find the neutron energy group given a neutron energy
unsigned int
ScatteringRecoilCrossSection::findNeutronEnergyGroup(Real energy)
{
  for (unsigned int g = 0; g < _G; ++g)
  {
    if (energy < _neutron_energy_limits[g] && energy > _neutron_energy_limits[g + 1])
      return g;
  }
  mooseError("Should never get here");
}

Real
ScatteringRecoilCrossSection::csvPrint(Real value) const
{
  if (std::abs(value) < _csv_tolerance)
    return 0;
  return value;
}

Real
ScatteringRecoilCrossSection::getSigmaRecoil(unsigned int g, unsigned int t, unsigned int l) const
{
  mooseAssert(g < _G, "g is larger than _G [indexed from g = 0, ..., G - 1]");
  mooseAssert(t < _T, "t is larger than _T [indexed from t = 0, ..., T - 1]");
  if (l <= _L)
    return _recoil_coefficient[l][t][g];
  return 0.0;
}

Real
ScatteringRecoilCrossSection::getMaxRecoilCosine(unsigned int g, unsigned t) const
{
  mooseAssert(g < _G, "g is larger than _G [indexed from g = 0, ..., G - 1]");
  mooseAssert(t < _T, "t is larger than _T [indexed from t = 0, ..., T - 1]");
  return _mu_L[t][g][1];
}

Real
ScatteringRecoilCrossSection::getMinRecoilCosine(unsigned int g, unsigned t) const
{
  mooseAssert(g < _G, "g is larger than _G [indexed from g = 0, ..., G - 1]");
  mooseAssert(t < _T, "t is larger than _T [indexed from t = 0, ..., T - 1]");
  return _mu_L[t][g][0];
}

bool
ScatteringRecoilCrossSection::isRecoilPossible(unsigned int g, unsigned int t) const
{
  return !(_mu_L[t][g][0] == -1 && _mu_L[t][g][1] == -1);
}

#endif
