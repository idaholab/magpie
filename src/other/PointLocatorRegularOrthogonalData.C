#include "PointLocatorRegularOrthogonalData.h"
#include "MooseError.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"

PointLocatorRegularOrthogonalData::PointLocatorRegularOrthogonalData(const std::vector<unsigned int> & cell_count,
                                                                     const Point & min_corner,
                                                                     const Point & max_corner,
                                                                     const MeshBase & mesh) :
    _cell_count(cell_count),
    _dim(_cell_count.size()),
    _min_corner(min_corner),
    _size(max_corner - min_corner)
{
  // determine total number of cells
  unsigned int num_root_elems = 1;
  for (auto count : _cell_count)
    num_root_elems *= count;

  // initialize cell storage
  _root_elems.assign(num_root_elems, nullptr);

  // sort all level 0 elements into the table
  auto  el = mesh.level_elements_begin(0);
  const auto end_el = mesh.level_elements_end(0);
  for (; el != end_el; ++el)
  {
    Elem * elem = *el;
    const unsigned int index = rootElementIndex(elem->centroid());
    _root_elems[index] = elem;
  }

  // sanity check
  for (auto el : _root_elems)
    if (el == nullptr)
      mooseError("Found a null root element");
}

const Elem *
PointLocatorRegularOrthogonalData::rootElement(const Point & p) const
{
  if (_root_elems[rootElementIndex(p)] == nullptr)
    mooseError("FCUK!");
  return _root_elems[rootElementIndex(p)];
}

unsigned int
PointLocatorRegularOrthogonalData::rootElementIndex(const Point & p) const
{
  unsigned int index = 0;
  for (unsigned int i = 0; i < _dim; ++i)
  {
    const int n = ((p(i) - _min_corner(i)) * _cell_count[i]) / _size(i);

    // check if we are inside the mesh TODO: respect _out_of_mesh_mode
    // mooseAssert(n >= 0 && n < int(_cell_count[i]), "Point is outside of the mesh");

    index = index * _cell_count[i] + (_cell_count[i] + n % _cell_count[i]) % _cell_count[i];
  }

  return index;
}
