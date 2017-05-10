#include "CoupledDefectAnnihilation.h"

template <>
InputParameters
validParams<CoupledDefectAnnihilation>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Coupled annihilation rate kernel that is proportional to the product of two coupled variables c and v.");
  params.addRequiredCoupledVar("c", "First defect concentration");
  params.addRequiredCoupledVar("v", "Other defect concentration");
  params.addParam<Real>("prefactor", 1.0, "Annihilation rate");
  return params;
}

CoupledDefectAnnihilation::CoupledDefectAnnihilation(const InputParameters & parameters)
  : Kernel(parameters),
    _c(coupledValue("c")),
    _c_var(coupled("c")),
    _v(coupledValue("v")),
    _v_var(coupled("v")),
    _prefactor(getParam<Real>("prefactor"))
{
}

Real
CoupledDefectAnnihilation::computeQpResidual()
{
  return _prefactor * _c[_qp] * _v[_qp] * _test[_i][_qp];
}

Real
CoupledDefectAnnihilation::computeQpJacobian()
{
  return 0.0;
}

Real
CoupledDefectAnnihilation::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_var)
    return _prefactor * _c[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  if (jvar == _c_var)
    return _prefactor * _phi[_j][_qp] * _v[_qp] * _test[_i][_qp];

  return 0.0;
}
