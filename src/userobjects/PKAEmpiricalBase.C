/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAEmpiricalBase.h"

template<>
InputParameters validParams<PKAEmpiricalBase>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  return params;
}

PKAEmpiricalBase::PKAEmpiricalBase(const InputParameters & parameters) :
    PKAGeneratorBase(parameters)
{
}

void
PKAEmpiricalBase::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real dt, Real vol, const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  mooseAssert(dt >= 0, "Passed a negative time window into PKAEmpiricalBase::appendPKAs");
  mooseAssert(vol >= 0, "Passed a negative volume into PKAEmpiricalBase::appendPKAs");

  const Real Z = getZ();
  const Real m = getM();
  const Real E = getE();

  int tag = ionTag(averaged_data._Z, averaged_data._M, Z, m);

  unsigned int num_pka = std::floor(dt * vol * getPKARate() + getRandomReal());

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

    // set location of the fission event
    setPosition(pka);

    // set random direction for ion 1 and opposite direction for ion 2
    setRandomDirection(pka);

    // add PKA to list
    ion_list.push_back(pka);
  }
}
