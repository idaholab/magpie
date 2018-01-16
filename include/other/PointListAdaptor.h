/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef POINTLISTADAPTOR_H
#define POINTLISTADAPTOR_H

#include "libmesh/nanoflann.hpp"
#include "libmesh/utility.h"

template <typename PointObject>
class PointListAdaptor
{
private:
  const std::vector<PointObject> & _pts;

public:
  PointListAdaptor(const std::vector<PointObject> & pts) : _pts(pts) {}

  /**
   * libMesh \p Point coordinate type
   */
  using coord_t = Real;

  /**
   * Must return the number of data points
   */
  inline size_t kdtree_get_point_count() const { return _pts.size(); }

  /**
   * get a Point reference from the PointData object at index idx in the list
   */
  const Point & getPoint(const size_t idx) const;

  /**
   * Returns the distance between the vector "p1[0:size-1]"
   * and the data point with index "idx_p2" stored in the class
   */
  inline coord_t kdtree_distance(const coord_t * p1, const size_t idx_p2, size_t /*size*/) const
  {
    mooseAssert(idx_p2 <= _pts.size(),
                "The point index should be less than"
                "total number of points used to build"
                "the KDTree.");

    const Point & p2 = getPoint(idx_p2);

    coord_t dist = 0.0;

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      dist += Utility::pow<2>(p1[i] - p2(i));

    return dist;
  }

  /**
   * Returns the dim'th component of the idx'th point in the class:
   * Since this is inlined and the "dim" argument is typically an immediate
   * value, the
   *  "if's" are actually solved at compile time.
   */
  inline coord_t kdtree_get_pt(const size_t idx, int dim) const
  {
    mooseAssert(dim < (int)LIBMESH_DIM,
                "The required component number should be less than the LIBMESH_DIM.");
    mooseAssert(idx < _pts.size(),
                "The index of the point should be less"
                "than total number of points used to"
                "construct the KDTree.");

    const Point & p = getPoint(idx);

    return p(dim);
  }

  /**
   * Optional bounding-box computation. This function is called by the nanoflann library.
   * If the return value is false, the standard bbox computation loop in the nanoflann
   * library is activated.
   */
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /* bb */) const
  {
    return false;
  }
};

// Specialization for PointListAdaptor<Point> (provide your own for custom types)
template <>
inline const Point &
PointListAdaptor<Point>::getPoint(const size_t idx) const { return _pts[idx]; }

#endif // POINTLISTADAPTOR_H
