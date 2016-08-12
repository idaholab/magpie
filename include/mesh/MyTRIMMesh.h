#ifndef MYTRIMMESH_H
#define MYTRIMMESH_H

#include "GeneratedMesh.h"
#include "PointLocatorRegularOrthogonal.h"

class MyTRIMMesh;

template<>
InputParameters validParams<MyTRIMMesh>();

class MyTRIMMesh : public GeneratedMesh
{
public:
  MyTRIMMesh(const InputParameters & parameters);
  MyTRIMMesh(const MyTRIMMesh & other_mesh);

  // No copy assignment
  MooseMesh & operator=(const MooseMesh & other_mesh) = delete;
  virtual MooseMesh & clone() const override;

  /// obtain a specialized PointLocator for the current mesh
  virtual std::unique_ptr<PointLocatorBase> getPointLocator() const override;

protected:
  std::vector<unsigned int> _cell_count;
  Point _min_corner;
  Point _max_corner;

  mutable std::unique_ptr<PointLocatorRegularOrthogonal> _point_locator;
};

#endif //MYTRIMMESH_H
