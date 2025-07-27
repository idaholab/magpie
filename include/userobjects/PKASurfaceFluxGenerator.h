/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "PKAGeneratorBase.h"

/**
 * Starts PKAs at a fixed point in a fixed direction
 */
class PKASurfaceFluxGenerator : public PKAGeneratorBase
{
public:
  static InputParameters validParams();

  PKASurfaceFluxGenerator(const InputParameters & parameters);

  virtual void appendPKAs(std::vector<MyTRIM_NS::IonBase> &,
                          const MyTRIMRasterizer::PKAParameters &,
                          const MyTRIMRasterizer::AveragedData &) const;

protected:
  /// provides a mean to override the angular distribution of the PKAs in derived class
  virtual void setDirection(MyTRIM_NS::IonBase & ion) const;

  /// the direction along which the PKAs start
  const RealVectorValue _direction;

  /// Uses point locator to determine the element id of the elemnt _point is in
  virtual void updateCachedElementID();

  /// number of PKAs to be generated per cm2-s
  const unsigned int _flux;

  /// time step
  const Real _dt;

  /// Mesh that comes from another generator
  const MeshBase * _mesh;

  /// boundary name
  const BoundaryName _boundary;

  /// surface area of boundary
  const Real _boundary_surface_area;

  /// the location from which to start PKAs
  Point _point;

  /// PKA nuclear charge
  const unsigned int _Z;

  /// PKA mass
  const Real _m;

  /// PKA Energy (in eV)
  const Real _E;

  /// cumulative probability list with element IDs
  const std::vector<std::pair<Real, dof_id_type>> _prob_elem_pairs;

  /// point locator to determine element pointers form locations
  std::unique_ptr<PointLocatorBase> _pl;

  /// the element id of the element containing _point
  dof_id_type _elem_id;

  Real boundarySurfaceArea(const BoundaryName & boundary, const std::unique_ptr<MeshBase> & mesh);

  std::vector<std::pair<Real, dof_id_type>>
  volumeWeightedElemDist(const BoundaryName & boundary, const std::unique_ptr<MeshBase> & mesh);
};
