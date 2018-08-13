/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAEmpiricalBase.h"

template <>
InputParameters
validParams<PKAEmpiricalBase>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  return params;
}

PKAEmpiricalBase::PKAEmpiricalBase(const InputParameters & parameters)
  : PKAGeneratorBase(parameters)
{
}

void
PKAEmpiricalBase::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                             const MyTRIMRasterizer::PKAParameters & pka_parameters,
                             const MyTRIMRasterizer::AveragedData &) const
{
  const auto dt = pka_parameters._dt;
  const auto vol = pka_parameters._volume;
  const auto recoil_rate_scaling = pka_parameters._recoil_rate_scaling;

  mooseAssert(dt >= 0, "Passed a negative time window into PKAEmpiricalBase::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAEmpiricalBase::appendPKAs");

  const Real Z = getZ();
  const Real m = getM();
  const Real E = getE();

  int tag = ionTag(pka_parameters, Z, m);

  unsigned int num_pka =
      std::floor(recoil_rate_scaling * dt * vol * getPKARate() + getRandomReal());

  for (unsigned i = 0; i < num_pka; ++i)
  {
    // each fission event generates a pair of recoils
    MyTRIM_NS::IonBase pka;

    // set charge, mass, energy
    pka._Z = Z;
    pka._m = m;
    pka._E = E;

    // the tag is the element this PKA get registered as upon stopping
    // -1 means the PKA will be ignored
    pka._tag = tag;

    // set stopping criteria
    pka.setEf();

    // set origin of the PKA
    setPosition(pka);

    // set a random direction for the PKA
    setRandomDirection(pka);

    // add PKA to list
    ion_list.push_back(pka);
  }
}
