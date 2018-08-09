/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef DEFECTANIHILATION_H
#define DEFECTANIHILATION_H

#include "Kernel.h"

class DefectAnnihilation;

template <>
InputParameters validParams<DefectAnnihilation>();

/**
 * Anihilation reaction kernel for vacancies and interstitials.
 */
class DefectAnnihilation : public Kernel
{
public:
  DefectAnnihilation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableValue & _v;
  unsigned int _v_var;
  const Real _prefactor;
};

#endif // DEFECTANIHILATION_H
