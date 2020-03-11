/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#include "ParkinCoulterInterface.h"
#include "PolyatomicDisplacementFunction.h"
#include "PolyatomicDamageEnergyFunction.h"
#include "PolyatomicDisplacementDerivativeFunction.h"
#include "MooseMesh.h"

// mytrim includes
#include <mytrim/element.h>

InputParameters
ParkinCoulterInterface::validParams()
{
  InputParameters params = emptyInputParameters();
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
  return params;
}

ParkinCoulterInterface::ParkinCoulterInterface(const MooseObject * moose_object)
  : _moose_obj(moose_object), _pars(moose_object->parameters())
{
  _Ecap = {{}};
  if (_pars.isParamValid("Ecap"))
    _Ecap = _pars.get<std::vector<std::vector<Real>>>("Ecap");
}

std::vector<MyTRIM_NS::Element>
ParkinCoulterInterface::polyMat() const
{
  std::vector<MyTRIM_NS::Element> poly_mat;
  std::vector<unsigned int> atomic_numbers = atomicNumbers();
  std::vector<Real> mass_numbers = massNumbers();
  std::vector<Real> N = numberFractions();
  std::vector<Real> threshold = _pars.get<std::vector<Real>>("displacement_thresholds");
  std::vector<Real> bind;
  if (_pars.isParamValid("lattice_binding_energies"))
    bind = _pars.get<std::vector<Real>>("lattice_binding_energies");
  else
    bind.assign(atomic_numbers.size(), 0.0);

  // perform some checks
  if (atomic_numbers.size() != mass_numbers.size() || atomic_numbers.size() != N.size() ||
      atomic_numbers.size() != threshold.size() || atomic_numbers.size() != bind.size())
    mooseError("Size mismatch for at least one parameter array. Z, A, number_fraction, "
               "displacement_thresholds and lattice_binding_energies"
               "must all have the same length.");

  for (unsigned int j = 0; j < atomic_numbers.size(); ++j)
  {
    MyTRIM_NS::Element element;
    element._Z = atomic_numbers[j];
    element._m = mass_numbers[j];
    element._t = N[j];
    element._Edisp = threshold[j];
    element._Elbind = bind[j];
    poly_mat.push_back(element);
  }
  return poly_mat;
}

void
ParkinCoulterInterface::computeDamageFunctions()
{
  // callback for allocating damage functions
  initDamageFunctions();

  Real energy = _padf->minEnergy();
  Real Emax = maxEnergy();
  Real threshold = _pars.get<Real>("uniform_energy_spacing_threshold");
  Real dE = _pars.get<Real>("uniform_energy_spacing");
  Real logdE = _pars.get<Real>("logarithmic_energy_spacing");

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

#endif
