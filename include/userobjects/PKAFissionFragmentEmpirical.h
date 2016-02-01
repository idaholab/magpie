#ifndef PKAFISSIONFRAGMENTEMPIRICAL_H
#define PKAFISSIONFRAGMENTEMPIRICAL_H

#include "PKAGeneratorBase.h"
#include "mytrim/invert.h"

class PKAFissionFragmentEmpirical;

template<>
InputParameters validParams<PKAFissionFragmentEmpirical>();

/**
 *
 */
class PKAFissionFragmentEmpirical : public PKAGeneratorBase
{
public:
  PKAFissionFragmentEmpirical(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::ionBase> & ion_list, Real dt, Real vol) const = 0;

protected:
  /// Fission rate (per unit volume) assuming pure fully dense UO2
  const Real _fission_rate;

  /**
   * Variable for the relative Uranium density (0..1).
   * Use a CONSTANT MONOMIAL AuxVariable for this.
   */
  const VariableValue & _relative_density;

  /// Mass inverter to sample PKA mass distribution
  MyTRIM_NS::massInverter _mass_inverter;

  /// Energy inverter to sample PKA energy distribution
  MyTRIM_NS::energyInverter _energy_inverter;
};

#endif // PKAFISSIONFRAGMENTEMPIRICAL_H
