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

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData &) const;

protected:
  /// Fission rate (per unit volume) assuming pure fully dense UO2
  const Real _fission_rate;

  /**
   * Variable for the relative Uranium density (0..1).
   * Use a CONSTANT MONOMIAL AuxVariable for this.
   */
  const VariableValue & _relative_density;
};

#endif // PKAFISSIONFRAGMENTEMPIRICAL_H
