#ifndef MYTRIMMESH_H
#define MYTRIMMESH_H

#include "GeneratedMesh.h"
#include "UnsignedIntegerPoint.h"

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

  virtual void buildMesh() override;

  ///@{ Public interface methods to get the mesh structure
  const UnsignedIntegerPoint & getCellCount() const { return _cell_count; }
  const Point & minCorner() const { return _min_corner; }
  const Point & maxCorner() const { return _max_corner; }
  ///@}

protected:
  UnsignedIntegerPoint _cell_count;
  Point _min_corner;
  Point _max_corner;
};

#endif //MYTRIMMESH_H
