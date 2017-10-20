/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAConstant.h"

template<>
InputParameters validParams<PKAConstant>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addParam<Real>("pka_rate", 1e-8, "PKA rate per unit volume (uses mesh units defined in the rasterizer and moose time units)");
  params.addRequiredParam<Real>("Z", "PKA nuclear charge");
  params.addRequiredParam<Real>("m", "PKA mass in amu");
  params.addRequiredParam<Real>("E", "PKA energy in eV");
  return params;
}

PKAConstant::PKAConstant(const InputParameters & parameters) :
    PKAGeneratorBase(parameters),
    _pka_rate(getParam<Real>("pka_rate")),
    _Z(getParam<Real>("Z")),
    _m(getParam<Real>("m")),
    _E(getParam<Real>("E"))
{
}

void
PKAConstant::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  mooseAssert(dt >= 0, "Passed a negative time window into PKAConstant::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAConstant::appendPKAs");

  int tag = ionTag(averaged_data._Z, averaged_data._M, _Z, _m);

  unsigned int num_pka = std::floor(dt * vol * _pka_rate + getRandomReal());

  for (unsigned i = 0; i < num_pka; ++i)
  {
    // each fission event generates a pair of recoils
    MyTRIM_NS::IonBase pka;

    // sample fission fragment masses
    pka._Z = _Z;
    pka._m = _m;
    pka._E = _E;

    // the tag is the element this PKA get registered as upon stopping
    // -1 means the PKA will be ignored
    pka._tag = tag;

    // set stopping criteria
    pka.setEf();

    // set location of the fission event
    setPosition(pka);

    // set random direction for ion 1 and opposite direction for ion 2
    setRandomDirection(pka);

    // add PKA to list
    ion_list.push_back(pka);
  }
}
