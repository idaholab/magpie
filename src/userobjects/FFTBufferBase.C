/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "FFTBufferBase.h"
#include "MyTRIMMesh.h"

#include "MooseTypes.h"
#include "AuxiliarySystem.h"
#include "AllLocalDofIndicesThread.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

#include <type_traits>

template <typename T>
InputParameters
FFTBufferBase<T>::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription("Fourier transform data buffer object.");
  params.addRangeCheckedParam<std::vector<int>>(
      "grid",
      "grid > 0",
      "Number of grid cells in each dimension to compute "
      "the FFT on (can be omitted when using MyTRIMMesh)");
  params.addCoupledVar("moose_variable",
                       "Optional AuxVariable to take the initial condition of the buffer from.");

  // make sure we run the object on initial to apply the initial conditions
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;

  return params;
}

template <typename T>
FFTBufferBase<T>::FFTBufferBase(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _dim(_mesh.dimension()),
    _cell_volume(1.0),
    _buffer_size(1),
    _moose_variable(coupledComponents("moose_variable")),
    _how_many(howMany())
{
  // make sure Real is double
  static_assert(
      std::is_same<Real, double>::value,
      "Libmesh must be build with double precision floating point numbers as the 'Real' type.");

  // FFT needs a regular grid of values to operate on
  if (isParamValid("grid"))
  {
    // cannot specify a corresponding MOOSE Variable when running mesh-less
    if (!_moose_variable.empty())
      paramError("moose_variable",
                 "You cannot specify a corresponding MOOSE Variable when running mesh-less, i.e. "
                 "using the 'grid' parameter.");

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

  // check optional coupled MOOSE variable (only for Real valued buffers)
  if (!_moose_variable.empty())
  {
    if (_moose_variable.size() != _how_many)
      paramError("moose_variable", "Variable needs to have ", _how_many, " components.");

    for (unsigned i = 0; i < _moose_variable.size(); ++i)
    {
      // check
      auto var = getVar("moose_variable", i);
      if (var->order() != 0)
        paramError("moose_variable", "Variable needs to have CONSTANT order.");
      if (var->isNodal())
        paramError("moose_variable", "Variable must be elemental.");

      _moose_variable[i] = &coupledValue("moose_variable");
    }
  }

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
  _start = reinterpret_cast<Real *>(start(0));
  std::ptrdiff_t istride = reinterpret_cast<char *>(start(1)) - reinterpret_cast<char *>(_start);
  if (istride % sizeof(Real) != 0)
    mooseError("Invalid data alignment");
  istride /= sizeof(Real);
}

template <>
void
FFTBufferBase<Real>::execute()
{
  std::cout << 'A';
  // get  grid / buffer location
  Point centroid = _current_elem->centroid();
  std::size_t a = 0;
  for (unsigned int i = 0; i < _dim; ++i)
    a = a * _grid[i] + std::floor(((centroid(i) - _min_corner(i)) * _grid[i]) / _box_size(i));

  // copy solution value from the variables into the buffer
  for (unsigned i = 0; i < _moose_variable.size(); ++i)
    _start[a * _how_many + i] = (*_moose_variable[i])[0];
}

template <typename T>
const T &
FFTBufferBase<T>::operator()(const Point & p) const
{
  std::size_t a = 0;
  for (unsigned int i = 0; i < _dim; ++i)
    a = a * _grid[i] + std::floor(((p(i) - _min_corner(i)) * _grid[i]) / _box_size(i));
  return _buffer[a];
}

template <typename T>
T &
FFTBufferBase<T>::operator()(const Point & p)
{
  std::size_t a = 0;
  for (unsigned int i = 0; i < _dim; ++i)
    a = a * _grid[i] + std::floor(((p(i) - _min_corner(i)) * _grid[i]) / _box_size(i));
  return _buffer[a];
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
