#ifndef COUPLEDDEFECTANIHILATION_H
#define COUPLEDDEFECTANIHILATION_H

#include "Kernel.h"

class CoupledDefectAnnihilation;

template <>
InputParameters validParams<CoupledDefectAnnihilation>();

/**
 * Coupled anihilation reaction kernel for vacancies and interstitials.
 */
class CoupledDefectAnnihilation : public Kernel
{
public:
  CoupledDefectAnnihilation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  ///@{ First reactant
  const VariableValue & _c;
  unsigned int _c_var;
  ///@}

  ///@{ Second reactant
  const VariableValue & _v;
  unsigned int _v_var;
  ///@}

  const Real _prefactor;
};

#endif // COUPLEDDEFECTANIHILATION_H
