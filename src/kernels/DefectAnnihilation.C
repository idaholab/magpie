/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "DefectAnnihilation.h"

template <>
InputParameters
validParams<DefectAnnihilation>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Annihilation rate kernel that is proportional to the product of the kernel variable and the coupled variable");
  params.addRequiredCoupledVar("v", "Other defect concentration");
  params.addParam<Real>("prefactor", 1.0, "Annihilation rate");
  return params;
}

DefectAnnihilation::DefectAnnihilation(const InputParameters & parameters)
  : Kernel(parameters),
    _v(coupledValue("v")),
    _v_var(coupled("v")),
    _prefactor(getParam<Real>("prefactor"))
{
}

Real
DefectAnnihilation::computeQpResidual()
{
  return _prefactor * _u[_qp] * _v[_qp] * _test[_i][_qp];
}

Real
DefectAnnihilation::computeQpJacobian()
{
  return _prefactor * _phi[_j][_qp] * _v[_qp] * _test[_i][_qp];
}

Real
DefectAnnihilation::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_var)
    return _prefactor * _u[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
