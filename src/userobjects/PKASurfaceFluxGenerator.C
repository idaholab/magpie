/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PKASurfaceFluxGenerator.h"
#include "MagpieUtils.h"

registerMooseObject("MagpieApp", PKASurfaceFluxGenerator);

InputParameters
PKASurfaceFluxGenerator::validParams()
{
  InputParameters params = PKAGeneratorBase::validParams();
  params.addClassDescription(
      "This PKAGenerator starts particles from a random point in a boundary in a fixed direction.");
  params.addRequiredParam<RealVectorValue>("direction", "The fixed direction the PKAs move along");
  params.addRequiredParam<BoundaryName>("boundary", "The boundary to apply the flux to.");
  params.addRequiredParam<unsigned int>(
      "flux",
      "The number of ions (starting PKAs) to strike the surface in ions/Area-second. Units of Area "
      "must be consistent with either mesh unit length or given boundary_surface_area units.");
  params.addRequiredParam<Real>("dt", "Time step size used in time stepper. Value in seconds.");
  params.addParam<Real>(
      "boundary_surface_area",
      "Surface area of boundary. Used to determine number of initial "
      "PKAs generated per timestep. Area must be consistent with flux dimensions.");
  params.addRequiredParam<Real>("Z", "PKA nuclear charge");
  params.addRequiredParam<Real>("m", "PKA mass in amu");
  params.addRequiredParam<Real>("E", "PKA energy in eV");
  return params;
}

PKASurfaceFluxGenerator::PKASurfaceFluxGenerator(const InputParameters & parameters)
  : PKAGeneratorBase(parameters),
    _direction(getParam<RealVectorValue>("direction")),
    _flux(getParam<unsigned int>("flux")),
    _dt(getParam<Real>("dt")),
    _mesh(_fe_problem.mesh()),
    _boundary(getParam<BoundaryName>("boundary")),
    _boundary_surface_area(PKASurfaceFluxGenerator::boundarySurfaceArea(_boundary)),
    _Z(getParam<Real>("Z")),
    _m(getParam<Real>("m")),
    _E(getParam<Real>("E")),
    _prob_elem_pairs(PKASurfaceFluxGenerator::volumeWeightedElemDist(_boundary)),
    _pl(_mesh.getPointLocator())
{
  _pl->enable_out_of_mesh_mode();
  updateCachedElementID();
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
    num_pka = std::floor(pka_parameters._recoil_rate_scaling * num_pka + getRandomReal());

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

    // sample elements
    Real rnd_num = getRandomReal();
    auto it = std::lower_bound(_prob_elem_pairs.begin(),
                               _prob_elem_pairs.end(),
                               rnd_num,
                               [](const std::pair<Real, dof_id_type> & entry, const Real & value)
                               { return entry.first < value; });

    // Get the selected element ID
    dof_id_type selected_elem_id = it->second;

    // Retrieve the element pointer from mesh
    const Elem * rnd_elem = _mesh.getMesh().elem_ptr(selected_elem_id);

    Point _point = MagpieUtils::randomElementPoint(
        *rnd_elem, getRandomPoint()); // sample random point in element
    pka._pos = _point;

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

Real
PKASurfaceFluxGenerator::boundarySurfaceArea(const BoundaryName & boundary)
{
  if (isParamValid("boundary_surface_area"))
    return getParam<Real>("boundary_surface_area");

  BoundaryID boundary_id = _mesh.getBoundaryID(boundary);
  Real volume_sum = 0;

  const auto range = _mesh.getBoundaryElementRange();
  for (const BndElement * bnd_elem : *range)
  {
    if (bnd_elem->_bnd_id != boundary_id)
      continue; // skip other boundaries
    const auto elem = bnd_elem->_elem;
    const auto side = bnd_elem->_side;
    volume_sum += elem->side_ptr(side)->volume();
  }
  mooseAssert(volume_sum > libMesh::TOLERANCE, "boundary_surface_area is not strictly positive!");
  return volume_sum;
}

std::vector<std::pair<Real, dof_id_type>>
PKASurfaceFluxGenerator::volumeWeightedElemDist(const BoundaryName & boundary)
{
  BoundaryID boundary_id = _mesh.getBoundaryID(boundary);
  Real volume_sum = 0;
  std::vector<std::pair<dof_id_type, Real>> elem_volumes; // (elem_id, side_volume)
  std::vector<std::pair<Real, dof_id_type>> prob_elem_pairs;

  const auto range = _mesh.getBoundaryElementRange();
  for (const BndElement * bnd_elem : *range)
  {
    if (bnd_elem->_bnd_id != boundary_id)
      continue; // skip other boundaries
    const auto elem = bnd_elem->_elem;
    dof_id_type elem_id = elem->id();
    const auto side = bnd_elem->_side;
    Real side_volume = elem->side_ptr(side)->volume();
    volume_sum += side_volume;
    elem_volumes.emplace_back(elem_id, side_volume);
  }

  // Normalize into cumulative distribution
  Real cumulative_sum = 0.0;
  for (const auto & [elem_id, side_vol] : elem_volumes)
  {
    cumulative_sum += side_vol / volume_sum;
    prob_elem_pairs.emplace_back(cumulative_sum, elem_id);
  }
  mooseAssert(std::abs(cumulative_sum - 1.0) < libMesh::TOLERANCE,
              "Volumes are not normalized to sum to 1.0!");
  mooseAssert(std::abs(prob_elem_pairs.back().first - 1.0) < libMesh::TOLERANCE,
              "Element probabilities are not normalized to sum to 1.0!");

  return prob_elem_pairs;
}

void
PKASurfaceFluxGenerator::updateCachedElementID()
{
  // get element containing the point
  mooseAssert(_pl != nullptr, "initialize() must be called on the MooseMyTRIMSample object.");
  const Elem * elem = (*_pl)(_point);
  if (elem == nullptr)
    mooseError("Point ", _point, " is not within the domain.");
  _elem_id = elem->id();
}
