/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PointLocatorRegularOrthogonalData.h"
#include "MooseError.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"

PointLocatorRegularOrthogonalData::PointLocatorRegularOrthogonalData(
    const std::vector<unsigned int> & cell_count,
    const Point & min_corner,
    const Point & max_corner,
    const MeshBase & mesh)
  : _dim(cell_count.size()),
    _cell_count(cell_count),
    _min_corner(min_corner),
    _cell_size(max_corner - min_corner)
{
  // determine total number of cells and their size
  unsigned int num_root_elems = 1;
  for (unsigned int i = 0; i < _dim; ++i)
  {
    if (_cell_count[i] == 0)
      mooseError("Cannot have zero cells in any spatial direction");

    num_root_elems *= _cell_count[i];
    _cell_size(i) /= _cell_count[i];
  }

  // initialize cell storage
  _root_elems.assign(num_root_elems, nullptr);

  // sort all level 0 elements into the table
  auto el = mesh.level_elements_begin(0);
  const auto end_el = mesh.level_elements_end(0);
  for (; el != end_el; ++el)
  {
    const Elem * elem = *el;
    const unsigned int index = rootElementIndex(elem->vertex_average());
    _root_elems[index] = elem;
  }

  // sanity check
  for (auto el : _root_elems)
    if (el == nullptr)
      mooseError("Found a null root element");
}

const Elem *
PointLocatorRegularOrthogonalData::rootElement(const Point & p, Point & el_pos) const
{
  unsigned int index = 0;
  const Point p0 = p - _min_corner;

  for (unsigned int i = 0; i < _dim; ++i)
  {
    const int n = p0(i) / _cell_size(i);

    // outside of mesh, return null
    if (n < 0 || n >= _cell_count[i])
      return nullptr;

    // cell internal coordinates from 0..1 used for bisecting refined cells
    el_pos(i) = p0(i) / _cell_size(i) - n;

    // construct root element index
    index = index * _cell_count[i] + n;
  }

  return _root_elems[index];
}

unsigned int
PointLocatorRegularOrthogonalData::rootElementIndex(const Point & p) const
{
  unsigned int index = 0;
  for (unsigned int i = 0; i < _dim; ++i)
  {
    const int n = (p(i) - _min_corner(i)) / _cell_size(i);

    if (n < 0 || n >= _cell_count[i])
      mooseError("Point not inside regular orthogonal mesh");

    index = index * _cell_count[i] + n;
  }

  return index;
}
