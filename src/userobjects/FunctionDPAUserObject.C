/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef GSL_ENABLED

#include "FunctionDPAUserObject.h"
#include "PolyatomicDisplacementFunction.h"

registerMooseObject("MagpieApp", FunctionDPAUserObject);

InputParameters
FunctionDPAUserObject::validParams()
{
  InputParameters params = DPAUserObjectBase::validParams();
  params.addRequiredParam<std::vector<std::vector<FunctionName>>>(
      "damage_functions", "Damage functions for each combinations of projectiles and targets.");
  params.addParam<Real>(
      "max_energy_step_size",
      100,
      "The maximum energy step size used for integration and interpolation. Default is 100 eV.");
  params.addClassDescription("Computes the dose in dpa from composition, cross section, damage "
                             "type, and neutron flux. The damage functions are provided as MOOSE "
                             "functions where energy is provided via the time arg slot.");
  return params;
}

FunctionDPAUserObject::FunctionDPAUserObject(const InputParameters & parameters)
  : DPAUserObjectBase(parameters), _max_delta_E(getParam<Real>("max_energy_step_size"))
{
  // get the damage functions
  std::vector<std::vector<FunctionName>> nm =
      getParam<std::vector<std::vector<FunctionName>>>("damage_functions");
  _damage_functions.resize(nm.size());
  for (unsigned int j = 0; j < nm.size(); ++j)
  {
    _damage_functions[j].resize(nm[j].size());
    for (unsigned int i = 0; i < nm[j].size(); ++i)
      _damage_functions[j][i] = &getFunctionByName(nm[j][i]);
  }
}

void
FunctionDPAUserObject::initialSetup()
{
  prepare();

  // make sure the dimension of _damage_functions is correct
  if (_damage_functions.size() != _atomic_numbers.size())
    paramError("damage_functions",
               "Size of leading dimension ",
               _damage_functions.size(),
               " is not identical to number of isotopes ",
               _atomic_numbers.size());
  for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
    if (_damage_functions[j].size() != _atomic_numbers.size())
      paramError("damage_functions",
                 "Size of entry ",
                 j,
                 " is ",
                 _damage_functions[j].size(),
                 " which is not identical to number of isotopes ",
                 _atomic_numbers.size());
}

void
FunctionDPAUserObject::execute()
{
  accumulateDamage();
}

void
FunctionDPAUserObject::onCompositionChanged()
{
  // compute the integral damage functions
  _integral_damage_functions.resize(_atomic_numbers.size());
  for (unsigned int i = 0; i < _atomic_numbers.size(); ++i)
  {
    // resize integral damage function
    _integral_damage_functions[i].resize(_atomic_numbers.size());

    for (unsigned int j = 0; j < _atomic_numbers.size(); ++j)
    {
      // loop over all energies
      Real energy = 0;
      Real max_energy = getMaxEnergy();
      std::vector<Real> energies(1);
      std::vector<Real> integral(1);
      // this is the energy time step, initially it must be
      // small but then can up to _max_delta_E
      Real deltaE = 1e-5;
      bool keep_going = true;
      while (keep_going)
      {
        Real energy_old = energy;
        energy += deltaE;
        if (energy >= max_energy)
        {
          energy = max_energy;
          keep_going = false;
        }

        // incremental integration from energy_old to energy
        std::vector<Real> points;
        std::vector<Real> weights;
        PolyatomicDisplacementFunction::gslQuadRule(100, points, weights);

        Real tint = 0;
        for (unsigned int qp = 0; qp < points.size(); ++qp)
        {
          Real e = 0.5 * (points[qp] + 1) * deltaE + energy_old;
          Real w = 0.5 * weights[qp] * deltaE;
          tint += w * _damage_functions[i][j]->value(e, Point());
        }

        // update energies and integral vectors
        energies.push_back(energy);
        Real old_integral = integral.back();
        integral.push_back(tint + old_integral);

        // now do an adjustment of the timestep but only up to _max_delta_E
        // the factor of 1.01 was found to be well converged for verification_2.i
        Real proposed_deltaE = deltaE * 1.01;
        if (proposed_deltaE < _max_delta_E)
          deltaE = proposed_deltaE;
      }

      // now the linear interpolation object is constructed and stored in integral damage function
      // var
      _integral_damage_functions[i][j] = LinearInterpolation(energies, integral);
    }
  }
}

Real
FunctionDPAUserObject::integralDamageFunction(Real T, unsigned int i, unsigned int j) const
{
  return _integral_damage_functions[i][j].sample(T);
}

void
FunctionDPAUserObject::finalize()
{
}

#endif
