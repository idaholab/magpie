/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPFMInterfaceIC.h"
#include "libmesh/utility.h"

registerMooseObject("MagpieApp", PolarPFMInterfaceIC);

InputParameters
PolarPFMInterfaceIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addClassDescription("Upsilon IC, eq. (50) from Physical Review B 89, 184102 (2014)");
  params.addRequiredCoupledVar("theta", "Theta order parameter");
  params.addRequiredParam<Real>("a_beta", "Interpolation coefficient a_beta");
  params.addRequiredParam<Real>("beta10", "Gradient energy coefficient between solid 1 and melt");
  params.addRequiredParam<Real>("beta20", "Gradient energy coefficient between solid 2 and melt");
  params.addParam<Real>("p", 1.0, "Interface parameter");
  return params;
}

PolarPFMInterfaceIC::PolarPFMInterfaceIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _theta(coupledValue("theta")),
    _a_beta(getParam<Real>("a_beta")),
    _beta10(getParam<Real>("beta10")),
    _beta20(getParam<Real>("beta20")),
    _p(getParam<Real>("p")),
    _xmin(_fe_problem.mesh().getMinInDimension(0)),
    _xmax(_fe_problem.mesh().getMaxInDimension(0))
{
}

Real
PolarPFMInterfaceIC::value(const Point & r)
{
  /// switching function
  const Real q = _a_beta * Utility::pow<2>(_theta[_qp]) -
                 2.0 * (_a_beta - 2.0) * Utility::pow<3>(_theta[_qp]) +
                 (_a_beta - 3.0) * Utility::pow<4>(_theta[_qp]);

  // solid-melt gradient energy coefficient (7)
  const Real betaS0 = _beta10 + (_beta20 - _beta10) * q;

  // sample width and position
  const Real W = _xmax - _xmin;
  const Real x = r(0) - _xmin;

  const Real ds0 = _p * sqrt(betaS0);
  return 1.0 - 1.0 / (1.0 + std::exp(-_p / ds0 * (x - W / 4.0))) +
         1.0 / (1.0 + std::exp(-_p / ds0 * (x - 3 * W / 4.0)));
}
