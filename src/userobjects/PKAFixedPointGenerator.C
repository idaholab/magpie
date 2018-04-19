/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKAFixedPointGenerator.h"
#include "MooseMesh.h"

registerMooseObject("MagpieApp", PKAFixedPointGenerator);

template<>
InputParameters validParams<PKAFixedPointGenerator>()
{
  InputParameters params = validParams<PKAGeneratorBase>();
  params.addClassDescription("This PKAGenerator starts particle from a fixed point in a random direction (isotropically).");
  params.addParam<unsigned int>("num_pkas", 1000, "The number of PKAs to be started from this position");
  params.addRequiredParam<Point>("point", "The point from which the PKAs are started");
  params.addRequiredParam<Real>("Z", "PKA nuclear charge");
  params.addRequiredParam<Real>("m", "PKA mass in amu");
  params.addRequiredParam<Real>("E", "PKA energy in eV");
  return params;
}

PKAFixedPointGenerator::PKAFixedPointGenerator(const InputParameters & parameters) :
    PKAGeneratorBase(parameters),
    _num_pka(getParam<unsigned int>("num_pkas")),
    _point(getParam<Point>("point")),
    _Z(getParam<Real>("Z")),
    _m(getParam<Real>("m")),
    _E(getParam<Real>("E")),
    _pl(_mesh.getPointLocator())
{
  _pl->enable_out_of_mesh_mode();
  updateCachedElementID();
}

void
PKAFixedPointGenerator::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list, Real /*dt*/, Real /*vol*/, Real recoil_rate_scaling, const MyTRIMRasterizer::AveragedData & averaged_data) const
{
  if (_current_elem->id() != _elem_id)
    return;

  unsigned int num_pka = _num_pka;
  if (recoil_rate_scaling != 1)
    num_pka = std::floor(recoil_rate_scaling * _num_pka + getRandomReal());

  for (unsigned i = 0; i < _num_pka; ++i)
  {
    // each fission event generates a pair of recoils
    MyTRIM_NS::IonBase pka;

    // sample fission fragment masses
    pka._Z = _Z;
    pka._m = _m;
    pka._E = _E;

    // the tag is the element this PKA get registered as upon stopping
    // -1 means the PKA will be ignored
    pka._tag = ionTag(averaged_data._Z, averaged_data._M, pka._Z, pka._m);;

    // set stopping criteria
    pka.setEf();

    // set initial location of the PKAs
    pka._pos = _point;

    // set random direction for ion 1 and opposite direction for ion 2
    setDirection(pka);

    // add PKA to list
    ion_list.push_back(pka);
  }
}

void
PKAFixedPointGenerator::setDirection(MyTRIM_NS::IonBase & ion) const
{
  // by default the angular distribution of PKAs is just uniform
  setRandomDirection(ion);
}

void
PKAFixedPointGenerator::updateCachedElementID()
{
  // get element containing the point
  mooseAssert(_pl != nullptr, "initialize() must be called on the MooseMyTRIMSample object.");
  const Elem * elem = (*_pl)(_point);
  if (elem == nullptr)
    mooseError("Point ", _point, " is not within the domain.");
  _elem_id = elem->id();
}
