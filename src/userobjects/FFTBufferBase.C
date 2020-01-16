/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#include "FFTBufferBase.h"
#include "MooseTypes.h"
#include "MyTRIMMesh.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

#include <type_traits>

// template <>
// InputParameters
// validParams<FFTBufferBase>()
// {
//   InputParameters params = validParams<GeneralUserObject>();
//   params.addClassDescription("");
//   return params;
// }

template <typename T>
InputParameters
FFTBufferBase<T>::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Fourier transform data buffer object.");
  params.addRangeCheckedParam<std::vector<int>>(
      "grid",
      "grid > 0",
      "Number of grid cells in each dimension to compute "
      "the FFT on (can be omitted when using MyTRIMMesh)");
  params.addParam<unsigned int>(
      "max_h_level", 0, "Further grid refinement to apply when using MyTRIMMesh with adaptivity");
  return params;
}

template <typename T>
FFTBufferBase<T>::FFTBufferBase(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension()),
    _cell_volume(1.0),
    _buffer_size(1)

{
  // make sure Real is double
  static_assert(
      std::is_same<Real, double>::value,
      "Libmesh must be build with double precision floating point numbers as the 'Real' type.");

  // FFT needs a regular grid of values to operate on
  if (isParamValid("grid"))
  {
    // if valid use the user-supplied sampling grid
    _grid = getParam<std::vector<int>>("grid");
    if (_grid.size() != _dim)
      paramError("grid", "Number of entries must match the mesh dimension ", _dim);
  }
  else
  {
    auto * mytrim_mesh = dynamic_cast<MyTRIMMesh *>(&_mesh);
    if (!mytrim_mesh)
      mooseError("Must specify 'grid' parameter if not using MyTRIMMesh");

    // querying the mesh to set grid
    for (unsigned int i = 0; i < _dim; ++i)
      _grid.push_back(mytrim_mesh->getCellCountInDimension(i));
  }

  // refine grid by max_h_level powers of two
  auto max_h_level = getParam<unsigned int>("max_h_level");
  for (auto & count : _grid)
    count = count << max_h_level;

  // get mesh extents and calculate space required and estimate spectrum bins
  for (unsigned int i = 0; i < _dim; ++i)
  {
    _min_corner(i) = _mesh.getMinInDimension(i);
    _max_corner(i) = _mesh.getMaxInDimension(i);
    _box_size(i) = _max_corner(i) - _min_corner(i);
    _cell_volume *= _box_size(i) / _grid[i];

    // last direction needs to be padded for in-place transforms
    _buffer_size *= (i == _dim - 1) ? ((_grid[i] >> 1) + 1) << 1 : _grid[i];
  }
  _buffer.resize(_buffer_size);

  // compute stride and start pointer
  auto istart = start(0);
  std::ptrdiff_t istride = reinterpret_cast<char *>(start(1)) - reinterpret_cast<char *>(istart);
  if (istride % sizeof(Real) != 0)
    mooseError("Invalid data alignment");
  istride /= sizeof(Real);
}

template <>
Real *
FFTBufferBase<Real>::start(std::size_t i)
{
  return &_buffer[i];
}

template <>
Real *
FFTBufferBase<RealVectorValue>::start(std::size_t i)
{
  return &_buffer[i](0);
}

template <>
Real *
FFTBufferBase<RankTwoTensor>::start(std::size_t i)
{
  return &_buffer[i](0, 0);
}

template <>
Real *
FFTBufferBase<RankThreeTensor>::start(std::size_t i)
{
  return &_buffer[i](0, 0, 0);
}

template <>
Real *
FFTBufferBase<RankFourTensor>::start(std::size_t i)
{
  return &_buffer[i](0, 0, 0, 0);
}

template <>
std::size_t
FFTBufferBase<Real>::howMany() const
{
  return 1;
}

template <>
std::size_t
FFTBufferBase<RealVectorValue>::howMany() const
{
  return LIBMESH_DIM;
}

template <>
std::size_t
FFTBufferBase<RankTwoTensor>::howMany() const
{
  return LIBMESH_DIM * LIBMESH_DIM;
}

template <>
std::size_t
FFTBufferBase<RankThreeTensor>::howMany() const
{
  return LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM;
}

template <>
std::size_t
FFTBufferBase<RankFourTensor>::howMany() const
{
  return LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM;
}

// explicit instantiation
template class FFTBufferBase<Real>;
template class FFTBufferBase<RealVectorValue>;
template class FFTBufferBase<RankTwoTensor>;
template class FFTBufferBase<RankThreeTensor>;
template class FFTBufferBase<RankFourTensor>;
