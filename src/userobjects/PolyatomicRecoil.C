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

InputParameters
PolyatomicRecoil::validParams()
{
  InputParameters params = ParkinCoulterBase::validParams();
  params.addRequiredParam<std::vector<unsigned int>>("Z", "Atomic numbers");
  params.addRequiredParam<std::vector<Real>>("A", "Mass numbers");
  params.addRequiredParam<std::vector<Real>>("number_fraction", "Number fractions");
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
  params.addRequiredRangeCheckedParam<Real>(
      "Emax", "Emax > 0", "Maximum desired energy to which displacement functions are computed");
  params.addClassDescription(
      "PolyatomicRecoil allows computation of total and net displacement functions,"
      "damage energy functions, and the derivative of the net displacement functions w.r.t. number "
      "fractions.");
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

PolyatomicRecoil::PolyatomicRecoil(const InputParameters & parameters)
  : ParkinCoulterBase(parameters)
{
}

std::vector<unsigned int>
PolyatomicRecoil::atomicNumbers() const
{
  return getParam<std::vector<unsigned int>>("Z");
}

std::vector<Real>
PolyatomicRecoil::massNumbers() const
{
  return getParam<std::vector<Real>>("A");
}

std::vector<Real>
PolyatomicRecoil::numberFractions() const
{
  return getParam<std::vector<Real>>("number_fraction");
}

void
PolyatomicRecoil::initDamageFunctions()
{
  // set the displacement function type
  nrt_type type = TOTAL;
  if (getParam<MooseEnum>("damage_type") == "NET")
    type = NET;
  else if (getParam<MooseEnum>("damage_type") == "ENERGY")
    type = ENERGY;
  else if (getParam<MooseEnum>("damage_type") == "NET_DERIVATIVE")
    type = NET_DERIVATIVE;

  if (type == ENERGY)
    _padf = libmesh_make_unique<PolyatomicDamageEnergyFunction>(polyMat(), type, _Ecap);
  else if (type == NET_DERIVATIVE)
  {
    _padf = libmesh_make_unique<PolyatomicDisplacementFunction>(polyMat(), NET, _Ecap);
    _padf_derivative = libmesh_make_unique<PolyatomicDisplacementDerivativeFunction>(
        polyMat(), type, dynamic_cast<PolyatomicDisplacementFunction *>(_padf.get()), _Ecap);
  }
  else
    _padf = libmesh_make_unique<PolyatomicDisplacementFunction>(polyMat(), type, _Ecap);
}

void
PolyatomicRecoil::execute()
{
  computeDamageFunctions();
}

void
PolyatomicRecoil::finalize()
{
  std::vector<unsigned int> atomic_numbers = atomicNumbers();
  std::vector<Real> mass_numbers = massNumbers();

  // displacement functions
  if (isParamValid("displacement_file_base"))
  {
    std::ofstream displacement_file;
    displacement_file.open(getParam<std::string>("displacement_file_base") + ".csv");

    PolyatomicDisplacementFunction * displacement_function =
        dynamic_cast<PolyatomicDisplacementFunction *>(_padf.get());
    PolyatomicDamageEnergyFunction * energy_function =
        dynamic_cast<PolyatomicDamageEnergyFunction *>(_padf.get());

    displacement_file << "energy (eV)";
    if (_padf_derivative)
      for (unsigned int projectile = 0; projectile < _padf_derivative->nSpecies(); ++projectile)
        for (unsigned int target = 0; target < _padf_derivative->nSpecies(); ++target)
          for (unsigned int derivative = 0; derivative < _padf_derivative->nSpecies(); ++derivative)
          {
            unsigned int projectile_zaid =
                1000 * atomic_numbers[projectile] + mass_numbers[projectile];
            unsigned int target_zaid = 1000 * atomic_numbers[target] + mass_numbers[target];
            unsigned int derivative_zaid =
                1000 * atomic_numbers[derivative] + mass_numbers[derivative];
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
                1000 * atomic_numbers[projectile] + mass_numbers[projectile];
            unsigned int target_zaid = 1000 * atomic_numbers[target] + mass_numbers[target];
            displacement_file << "," << projectile_zaid << "->" << target_zaid;
          }
      else if (energy_function)
        for (unsigned int projectile = 0; projectile < _padf->nSpecies(); ++projectile)
        {
          unsigned int projectile_zaid =
              1000 * atomic_numbers[projectile] + mass_numbers[projectile];
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
                                       _padf_derivative->energyPoint(n),
                                       projectile,
                                       target,
                                       derivative);
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
                                                  _padf->energyPoint(n), projectile, target);
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
                              << energy_function->linearInterpolation(_padf->energyPoint(n),
                                                                      projectile);
          displacement_file << std::endl;
        }
      }
    }

    displacement_file.close();
  }
}

#endif
