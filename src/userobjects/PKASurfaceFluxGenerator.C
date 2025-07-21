/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKASurfaceFluxGenerator.h"

registerMooseObject("MagpieApp", PKASurfaceFluxGenerator);

InputParameters
PKASurfaceFluxGenerator::validParams()
{
  InputParameters params = PKAGeneratorBase::validParams();
  params.addClassDescription(
      "This PKAGenerator starts particles from a random point in a boundary in a fixed direction.");
  params.addRequiredParam<RealVectorValue>("direction", "The fixed direction the PKAs move along");
  params.addRequiredParam<BoundaryName>("boundary", "The boundary to apply the flux to.")
      params.addParam<unsigned int>(
          "flux", 1e8, "The number of ions (starting PKAs) to strike the surface in ions/cm2-sec");
  params.addRequiredParam<Real>("dt", "Time step size used in time stepper. Value in cm2!");
  params.addRequiredParam<Real>("boundary_surface_area",
                                "Surface area of boundary. Used to determine number of initial "
                                "PKAs generated per timestep. Value in cm2!");
  params.addRequiredParam<Real>("Z", "PKA nuclear charge");
  params.addRequiredParam<Real>("m", "PKA mass in amu");
  params.addRequiredParam<Real>("E", "PKA energy in eV");
  return params;
}

PKASurfaceFluxGenerator::PKASurfaceFluxGenerator(const InputParameters & parameters)
  : PKAFixedPointGenerator(parameters),
    _direction(getParam<RealVectorValue>("direction")),
    _flux(getParam<Real>("flux")),
    _boundary_surface_area<Real>(getParam("boundary_surface_area")),
    _dt<Real>(getParam("dt")),
    _boundary(getParam<BoundaryName>("boundary")),
    _Z(getParam<Real>("Z")),
    _m(getParam<Real>("m")),
    _E(getParam<Real>("E"))
{
}

void
PKASurfaceFluxGenerator::appendPKAs(std::vector<MyTRIM_NS::IonBase> & ion_list,
                                    const MyTRIMRasterizer::PKAParameters & pka_parameters,
                                    const MyTRIMRasterizer::AveragedData &) const
{
  if (_current_elem->id() != _elem_id)
    return;

  unsigned int num_pka = _flux * _dt * _boundary_surface_area;
  if (pka_parameters._recoil_rate_scaling != 1)
    num_pka = std::floor(pka_parameters._recoil_rate_scaling * _num_pka + getRandomReal());

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
    pka._tag = ionTag(pka_parameters, pka._Z, pka._m);
    ;

    // set stopping criteria
    pka.setEf();

    // set initial location of the PKAs
    pka._pos = (0, 0, 0); // CHANGE THIS

    // set random direction for ion 1 and opposite direction for ion 2
    setDirection(pka);

    // add PKA to list
    ion_list.push_back(pka);
  }
}

void
PKASurfaceFluxGenerator::setDirection(MyTRIM_NS::IonBase & ion) const
{
  ion._dir = _direction;
}
