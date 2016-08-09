#ifndef POINTLOCATORREGULARORTHOGONALDATA_H
#define POINTLOCATORREGULARORTHOGONALDATA_H

#include "Moose.h"
#include "libmesh/point.h"
#include "libmesh/reference_counted_object.h"

class PointLocatorRegularOrthogonalData;

namespace libMesh {
  class MeshBase;
  class Elem;
}

/**
 * Lookup data table for PointLocatorRegularOrthogonalData. This is only built
 * once for the master PointLocator that is held by the mesh. The mesh only hands out
 * sub locators that point to this structure on the master.
 * The lookup data is a regular orthogonal matrix of element pointers that correspond to
 * the level 0 elements in the regular orthogonal mesh.
 */
class PointLocatorRegularOrthogonalData :
  public ReferenceCountedObject<PointLocatorRegularOrthogonalData>
{
public:
  PointLocatorRegularOrthogonalData(const std::vector<unsigned int> & cell_count,
                                    const Point & min_corner,
                                    const Point & max_corner,
                                    const MeshBase & mesh);

  /// get the root element pointer
  const Elem * rootElement(const Point & p) const;

protected:
  /// get the root element index for the current point
  unsigned int rootElementIndex(const Point & p) const;

  /// number of cells (bins) along each dimension
  std::vector<unsigned int> _cell_count;

  /// underlying mesh dimension
  unsigned int _dim;

  /// mesh point with the smallest coordinate components (bottom left back corner)
  Point _min_corner;

  /// lengths of the sides of mesh
  Point _size;

  /// serialized nz*ny*nx list of pointers to the coarsest elements in the mesh
  std::vector<Elem *> _root_elems;
};

#endif //POINTLOCATORREGULARORTHOGONALDATA_H
