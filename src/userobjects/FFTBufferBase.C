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
    _moose_variable(coupledComponents("moose_variable")),
    _k_table(_dim),
    _how_many(_real_space_data.howMany())
{
  // make sure Real is double
  static_assert(
      std::is_same<Real, double>::value,
      "Libmesh must be build with double precision floating point numbers as the 'Real' type.");

  static_assert(sizeof(std::complex<Real>) == 2 * sizeof(Real),
                "Complex numbers based on 'Real' should have the size of two 'Real' types.");

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

      _moose_variable[i] = &coupledValue("moose_variable", i);
    }
  }

  // get mesh extents and calculate space required and estimate spectrum bins
  std::size_t real_space_buffer_size = 1;
  std::size_t reciprocal_space_buffer_size = 1;
  for (unsigned int i = 0; i < _dim; ++i)
  {
    _min_corner(i) = _mesh.getMinInDimension(i);
    _max_corner(i) = _mesh.getMaxInDimension(i);
    _box_size(i) = _max_corner(i) - _min_corner(i);
    _cell_volume *= _box_size(i) / _grid[i];

    // unumber of physical grid cells
    real_space_buffer_size *= _grid[i];
    // last direction needs to be roughly cut in half
    reciprocal_space_buffer_size *= (i == _dim - 1) ? ((_grid[i] >> 1) + 1) : _grid[i];

    // precompute kvector components for current direction
    _k_table[i].resize(_grid[i]);
    for (int j = 0; j < _grid[i]; ++j)
      _k_table[i][j] = 2.0 * libMesh::pi *
                       ((j < (_grid[i] >> 1) + 1) ? Real(j) : (Real(j) - _grid[i])) / _box_size(i);
  }
  _real_space_data.resize(real_space_buffer_size);
  _reciprocal_space_data.resize(reciprocal_space_buffer_size);

  // compute stride and start pointer
  _real_space_data_start = reinterpret_cast<Real *>(_real_space_data.start(0));
  _reciprocal_space_data_start = reinterpret_cast<Complex *>(_reciprocal_space_data.start(0));

  _real_space_data_stride = reinterpret_cast<char *>(_real_space_data.start(1)) -
                            reinterpret_cast<char *>(_real_space_data_start);
  _reciprocal_space_data_stride = reinterpret_cast<char *>(_reciprocal_space_data.start(1)) -
                                  reinterpret_cast<char *>(_reciprocal_space_data_start);

  if (_real_space_data_stride % sizeof(Real) != 0)
    mooseError("Invalid data alignment");
  _real_space_data_stride /= sizeof(Real);

  if (_reciprocal_space_data_stride % sizeof(Complex) != 0)
    mooseError("Invalid data alignment");
  _reciprocal_space_data_stride /= sizeof(Complex);
}

template <typename T>
void
FFTBufferBase<T>::execute()
{
  // get  grid / buffer location
  Point centroid = _current_elem->centroid();
  std::size_t a = 0;
  for (unsigned int i = 0; i < _dim; ++i)
    a = a * _grid[i] + std::floor(((centroid(i) - _min_corner(i)) * _grid[i]) / _box_size(i));

  // copy solution value from the variables into the buffer
  for (unsigned i = 0; i < _moose_variable.size(); ++i)
    _real_space_data_start[a * _how_many + i] = (*_moose_variable[i])[0];
}

template <typename T>
void
FFTBufferBase<T>::forward()
{
  forwardRaw();
}

template <typename T>
void
FFTBufferBase<T>::backward()
{
  backwardRaw();
  _real_space_data *= backwardScale();
}

template <typename T>
const T &
FFTBufferBase<T>::operator()(const Point & p) const
{
  std::size_t a = 0;
  for (unsigned int i = 0; i < _dim; ++i)
    a = a * _grid[i] + std::floor(((p(i) - _min_corner(i)) * _grid[i]) / _box_size(i));
  return _real_space_data[a];
}

template <typename T>
T &
FFTBufferBase<T>::operator()(const Point & p)
{
  std::size_t a = 0;
  for (unsigned int i = 0; i < _dim; ++i)
    a = a * _grid[i] + std::floor(((p(i) - _min_corner(i)) * _grid[i]) / _box_size(i));
  return _real_space_data[a];
}

// explicit instantiations
template class FFTBufferBase<Real>;
template class FFTBufferBase<RealVectorValue>;
template class FFTBufferBase<RankTwoTensor>;
template class FFTBufferBase<RankThreeTensor>;
template class FFTBufferBase<RankFourTensor>;
