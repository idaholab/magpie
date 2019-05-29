/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "PolyatomicRecoil.h"
#include "PolyatomicDisplacementFunction.h"
#include "PolyatomicDamageEnergyFunction.h"
#include "PolyatomicDisplacementDerivativeFunction.h"
#include "MooseMesh.h"

// mytrim includes
#include <mytrim/element.h>

registerMooseObject("MagpieApp", PolyatomicRecoil);

template <>
InputParameters
validParams<PolyatomicRecoil>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<std::vector<unsigned int>>("Z", "Atomic numbers");
  params.addRequiredParam<std::vector<Real>>("A", "Mass numbers");
  params.addRequiredParam<std::vector<Real>>("number_fraction", "Number fractions");
  params.addRequiredParam<std::vector<Real>>("displacement_thresholds", "Dispacement thresholds");
  params.addParam<std::vector<Real>>("lattice_binding_energies", "Lattice binding energies");
  params.addParam<std::vector<std::vector<Real>>>(
      "Ecap", "Capture energy Ecap_ij of species i being trapped in j site");
  params.addRangeCheckedParam<Real>("uniform_energy_spacing_threshold",
                                    10,
                                    "uniform_energy_spacing_threshold >= 0",
                                    "Threshold below which energy points are spaced uniformly.");
  params.addRangeCheckedParam<Real>("uniform_energy_spacing",
                                    0.25,
                                    "uniform_energy_spacing > 0",
                                    "Uniform energy spacing below the threshold");
  params.addRequiredRangeCheckedParam<Real>(
      "logarithmic_energy_spacing",
      "logarithmic_energy_spacing > 1",
      "Spacing of the energy points En in log space energy_spacing = E_{n+1} / En");
  params.addRequiredRangeCheckedParam<Real>(
      "Emax", "Emax > 0", "Maximum desired energy to which displacement functions are computed");
  MooseEnum nrt_damage_types("TOTAL NET ENERGY NET_DERIVATIVE", "TOTAL");
  params.addParam<MooseEnum>(
      "damage_type",
      nrt_damage_types,
      "NRT damage types TOTAL: total number of atoms that have been dispaced during cascade [nij "
      "in PK JNM, 101, 1981]\n"
      "NET: number of atoms displaced and not recaptured [nij in PK JNM, 101, 1981]\n"
      "ENERGY: energy deposted by recoil of type i and energy E[nu_i in PK JNM, 88, 1980]\n"
      "NET_DERIVATIVE: derivative of NET w.r.t. the partial number fractions.");
  params.addParam<std::string>("displacement_file_base",
                               "The output file base for displacement function in csv format.");
  params.addClassDescription(
      "PolyatomicRecoil allows computation of total and net displacement functions,"
      "damage energy functions, and the derivative of the net displacement functions w.r.t. number "
      "fractions.");
  return params;
}

PolyatomicRecoil::PolyatomicRecoil(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _atomic_numbers(getParam<std::vector<unsigned int>>("Z")),
    _mass_numbers(getParam<std::vector<Real>>("A"))
{
  std::vector<Real> N = getParam<std::vector<Real>>("number_fraction");
  std::vector<Real> threshold = getParam<std::vector<Real>>("displacement_thresholds");
  std::vector<Real> bind;
  if (isParamValid("lattice_binding_energies"))
    bind = getParam<std::vector<Real>>("lattice_binding_energies");
  else
    bind.assign(_atomic_numbers.size(), 0.0);

  std::vector<std::vector<Real>> Ecap = {{}};
  if (isParamValid("Ecap"))
    Ecap = getParam<std::vector<std::vector<Real>>>("Ecap");

  if (_atomic_numbers.size() != _mass_numbers.size() || _atomic_numbers.size() != N.size() ||
      _atomic_numbers.size() != threshold.size() || _atomic_numbers.size() != bind.size())
    mooseError("Size mismatch for at least one parameter array. Z, A, number_fraction, "
               "displacement_thresholds and lattice_binding_energies"
               "must all have the same length.");

  std::vector<MyTRIM_NS::Element> poly_mat;
  for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
  {
    MyTRIM_NS::Element element;
    element._Z = _atomic_numbers[j];
    element._m = _mass_numbers[j];
    element._t = N[j];
    element._Edisp = threshold[j];
    element._Elbind = bind[j];
    poly_mat.push_back(element);
  }

  // set the displacement function type
  nrt_type type = TOTAL;
  if (getParam<MooseEnum>("damage_type") == "NET")
    type = NET;
  else if (getParam<MooseEnum>("damage_type") == "ENERGY")
    type = ENERGY;
  else if (getParam<MooseEnum>("damage_type") == "NET_DERIVATIVE")
    type = NET_DERIVATIVE;

  if (type == ENERGY)
    _padf = libmesh_make_unique<PolyatomicDamageEnergyFunction>(poly_mat, type, Ecap);
  else if (type == NET_DERIVATIVE)
  {
    _padf = libmesh_make_unique<PolyatomicDisplacementFunction>(poly_mat, NET, Ecap);
    _padf_derivative = libmesh_make_unique<PolyatomicDisplacementDerivativeFunction>(
        poly_mat, type, dynamic_cast<PolyatomicDisplacementFunction *>(_padf.get()), Ecap);
  }
  else
    _padf = libmesh_make_unique<PolyatomicDisplacementFunction>(poly_mat, type, Ecap);

  if (getParam<Real>("Emax") < getParam<Real>("uniform_energy_spacing_threshold"))
    mooseError("Emax must be larger than uniform_energy_spacing_threshold.");
}

void
PolyatomicRecoil::execute()
{
  Real energy = _padf->minEnergy();
  Real Emax = getParam<Real>("Emax");
  Real threshold = getParam<Real>("uniform_energy_spacing_threshold");
  Real dE = getParam<Real>("uniform_energy_spacing");
  Real logdE = getParam<Real>("logarithmic_energy_spacing");

  for (;;) // while (energy <= Emax)
  {
    energy = energy < threshold ? energy + dE : energy * logdE;
    if (energy > Emax)
    {
      _padf->advanceDisplacements(Emax);
      break;
    }

    // increment displacements for value of energy
    _padf->advanceDisplacements(energy);
  }

  if (_padf_derivative)
    for (unsigned int n = 1; n < _padf->nEnergySteps(); ++n)
      _padf_derivative->advanceDisplacements(_padf->energyPoint(n));
}

void
PolyatomicRecoil::finalize()
{
  // displacement functions
  if (isParamValid("displacement_file_base"))
  {
    std::ofstream displacement_file;
    displacement_file.open(getParam<std::string>("displacement_file_base") + ".csv");

    const PolyatomicDisplacementFunction * displacement_function =
        dynamic_cast<PolyatomicDisplacementFunction *>(_padf.get());
    const PolyatomicDamageEnergyFunction * energy_function =
        dynamic_cast<PolyatomicDamageEnergyFunction *>(_padf.get());

    displacement_file << "energy (eV)";
    if (_padf_derivative)
      for (unsigned int projectile = 0; projectile < _padf_derivative->nSpecies(); ++projectile)
        for (unsigned int target = 0; target < _padf_derivative->nSpecies(); ++target)
          for (unsigned int derivative = 0; derivative < _padf_derivative->nSpecies(); ++derivative)
          {
            unsigned int projectile_zaid =
                1000 * _atomic_numbers[projectile] + _mass_numbers[projectile];
            unsigned int target_zaid = 1000 * _atomic_numbers[target] + _mass_numbers[target];
            unsigned int derivative_zaid =
                1000 * _atomic_numbers[derivative] + _mass_numbers[derivative];
            displacement_file << ", d(" << projectile_zaid << "->" << target_zaid << ") / d("
                              << derivative_zaid << ")";
          }
    else
    {
      if (displacement_function)
        for (unsigned int projectile = 0; projectile < _padf->nSpecies(); ++projectile)
          for (unsigned int target = 0; target < _padf->nSpecies(); ++target)
          {
            unsigned int projectile_zaid =
                1000 * _atomic_numbers[projectile] + _mass_numbers[projectile];
            unsigned int target_zaid = 1000 * _atomic_numbers[target] + _mass_numbers[target];
            displacement_file << "," << projectile_zaid << "->" << target_zaid;
          }
      else if (energy_function)
        for (unsigned int projectile = 0; projectile < _padf->nSpecies(); ++projectile)
        {
          unsigned int projectile_zaid =
              1000 * _atomic_numbers[projectile] + _mass_numbers[projectile];
          displacement_file << "," << projectile_zaid;
        }
      else
        mooseError("Displacement object does not have right type");
    }

    displacement_file << std::endl;

    if (_padf_derivative)
      for (unsigned int n = 0; n < _padf->nEnergySteps(); ++n)
      {
        displacement_file << _padf_derivative->energyPoint(n);
        for (unsigned int projectile = 0; projectile < _padf_derivative->nSpecies(); ++projectile)
          for (unsigned int target = 0; target < _padf_derivative->nSpecies(); ++target)
            for (unsigned int derivative = 0; derivative < _padf_derivative->nSpecies();
                 ++derivative)
              displacement_file << ","
                                << _padf_derivative->linearInterpolation(
                                       projectile,
                                       target,
                                       derivative,
                                       _padf_derivative->energyPoint(n));
        displacement_file << std::endl;
      }
    else
    {
      if (displacement_function)
      {
        for (unsigned int n = 0; n < _padf->nEnergySteps(); ++n)
        {
          displacement_file << _padf->energyPoint(n);
          for (unsigned int projectile = 0; projectile < _padf->nSpecies(); ++projectile)
            for (unsigned int target = 0; target < _padf->nSpecies(); ++target)
            {
              Real delta_ij =
                  getParam<MooseEnum>("damage_type") == "TOTAL" && target == projectile ? 1 : 0;
              displacement_file << ","
                                << delta_ij + displacement_function->linearInterpolation(
                                                  projectile, target, _padf->energyPoint(n));
            }
          displacement_file << std::endl;
        }
      }
      else if (energy_function)
      {
        for (unsigned int n = 0; n < _padf->nEnergySteps(); ++n)
        {
          displacement_file << _padf->energyPoint(n);
          for (unsigned int projectile = 0; projectile < _padf->nSpecies(); ++projectile)
            displacement_file << ","
                              << energy_function->linearInterpolation(projectile,
                                                                      _padf->energyPoint(n));
          displacement_file << std::endl;
        }
      }
    }

    displacement_file.close();
  }
}

#endif
