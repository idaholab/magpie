/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMMESH_H
#define MYTRIMMESH_H

#include "GeneratedMesh.h"
#include "PointLocatorRegularOrthogonal.h"

class MyTRIMMesh;

template <>
InputParameters validParams<MyTRIMMesh>();

/**
 * Restricted regular orthogonal generated mesh with equal size level 0 elements
 * (no bias) and EDGE2, QUAD4, or HEX8 elements only. These restrictions allow
 * the use of a custom super fast point locator to find elements based on spatial
 * locations (required for fast material property lookup in the Binary Collision
 * Monte Carlo stage). This mesh is also recommended for use with the FourierTransform
 * user object.
 */
class MyTRIMMesh : public GeneratedMesh
{
public:
  MyTRIMMesh(const InputParameters & parameters);
  MyTRIMMesh(const MyTRIMMesh & other_mesh);

  // No copy assignment
  MooseMesh & operator=(const MooseMesh & other_mesh) = delete;
  MooseMesh & clone() const override;

  /// obtain a specialized PointLocator for the current mesh
  std::unique_ptr<PointLocatorBase> getPointLocator() const override;

  Real getMinInDimension(unsigned int component) const override;
  Real getMaxInDimension(unsigned int component) const override;

  /// since this is a regular mesh we can report the grid size in each dimension
  unsigned int getCellCountInDimension(unsigned int component);

protected:
  std::vector<unsigned int> _cell_count;
  Point _min_corner;
  Point _max_corner;

  mutable std::unique_ptr<PointLocatorRegularOrthogonal> _point_locator;
};

#endif // MYTRIMMESH_H
