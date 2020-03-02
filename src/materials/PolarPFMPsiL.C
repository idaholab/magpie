/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPFMPsiL.h"

registerMooseObject("MagpieApp", PolarPFMPsiL);

InputParameters
PolarPFMPsiL::validParams()
{
  InputParameters params = DerivativeParsedMaterialHelper::validParams();
  params.addClassDescription("");
  params.addRequiredCoupledVar("upsilon", "Upsilon order parameter");
  params.addRequiredCoupledVar("theta", "Theta order parameter");
  params.addRequiredParam<Real>("G0", "Thermal energy of the melt");
  params.addRequiredParam<Real>("a_theta", "Interpolation coefficient a_theta");
  params.addRequiredParam<Real>("a_A", "Interpolation coefficient a_A");
  params.addRequiredParam<Real>("A10", "Barrier coefficient solid 1 and melt");
  params.addRequiredParam<Real>("A20", "Barrier coefficient solid 2 and melt");
  params.addRequiredParam<Real>("A21", "Barrier coefficient solid 2 and solid 1");
  params.addRequiredParam<Real>("DeltaG10",
                                "Difference in thermal energy between solid 1 and melt");
  params.addRequiredParam<Real>("DeltaG21",
                                "Difference in thermal energy between solid 2 and solid 1");
  return params;
}

PolarPFMPsiL::PolarPFMPsiL(const InputParameters & parameters)
  : DerivativeParsedMaterialHelper(parameters),
    _upsilon("upsilon"),
    _theta("theta"),
    _g0(getParam<Real>("G0")),
    _a_theta(getParam<Real>("a_theta")),
    _a_A(getParam<Real>("a_A")),
    _deltaG10(getParam<Real>("DeltaG10")),
    _deltaG21(getParam<Real>("DeltaG21")),
    _a10(getParam<Real>("A10")),
    _a20(getParam<Real>("A20")),
    _a21(getParam<Real>("A21"))
{
  // All equation numbers from Physical Review B 89, 184102 (2014)

  // interpolating function (8)
  EBFunction q;
  EBTerm y, a;
  q(y, a) = a * y * y - 2 * (a - 2) * y * y * y + (a - 3) * y * y * y * y;

  // change in thermal energy of the phases (5)
  EBTerm DeltaG = _deltaG10 + _deltaG21 * q(_theta, 0);

  // solid-melt energy barrier coefficient (6)
  EBTerm AS0 = _a10 + (_a20 - _a10) * q(_theta, _a_theta);

  // energy barrier (3)
  EBTerm psiB = AS0 * _upsilon * _upsilon * (1 - _upsilon) * (1 - _upsilon) +
                _a21 * _theta * _theta * (1 - _theta) * (1 - _theta) * q(_upsilon, _a_A);

  // thermal energy (2)
  EBTerm psi = _g0 + DeltaG * q(_upsilon, 0);

  functionParse(psi + psiB);
}
