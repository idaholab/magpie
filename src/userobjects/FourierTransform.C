/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#include "FourierTransform.h"
#include "MyTRIMMesh.h"

#include "libmesh/utility.h"

#include <functional>

registerMooseObject("MagpieApp", FourierTransform);

template <>
InputParameters
validParams<FourierTransform>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Compute the Fourier transform of a given variable field.");
  params.addCoupledVar("variable", "Variable field to compute the transform of");
  params.addRangeCheckedParam<std::vector<int>>(
      "grid",
      "grid > 0",
      "Number of grid cells in each dimension to compute "
      "the FFT on (can be omitted when using MyTRIMMesh)");
  params.addParam<unsigned int>(
      "max_h_level", 0, "Further grid refinement to apply when using MyTRIMMesh with adaptivity");
  return params;
}

FourierTransform::FourierTransform(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _var(coupledValue("variable")),
    _dim(_mesh.dimension()),
    _cell_volume(1.0),
    _buffer_size(1),
    _perf_plan(registerTimedSection("fftw_plan_r2r", 2)),
    _perf_fft(registerTimedSection("fftw_execute", 2))
{
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

  // initialize FFTW buffer and plan
  _buffer.resize(_buffer_size);
  std::vector<fftw_r2r_kind> kind(_dim, FFTW_R2HC);
  {
    TIME_SECTION(_perf_plan);
    _plan = fftw_plan_r2r(
        _dim, _grid.data(), _buffer.data(), _buffer.data(), kind.data(), FFTW_ESTIMATE);
  }
}

FourierTransform::~FourierTransform()
{
  // destroy FFTW plan
  fftw_destroy_plan(_plan);
}

void
FourierTransform::initialize()
{
  // avoid assign here to make _sure_ the data address (which has already been passed to the plan)
  // stays the same (i.e. no reallocations)!
  std::fill(_buffer.begin(), _buffer.end(), 0.0);
}

void
FourierTransform::execute()
{
  // loop over all quadrature points in the element
  for (std::size_t qp = 0; qp < _qrule->n_points(); ++qp)
  {
    // get FFT grid cell index for current QP
    std::size_t index = 0;
    for (unsigned int component = 0; component < _dim; ++component)
    {
      const unsigned int cell =
          ((_q_point[qp](component) - _min_corner(component)) * _grid[component]) /
          _box_size(component);
      index = index * _grid[component] + cell;
    }

    // add qp contribution to FFT grid cell
    _buffer[index] += _var[qp] * _JxW[qp] * _coord[qp] / _cell_volume;
  }
}

void
FourierTransform::finalize()
{
  // get the sum of all buffers to all processes
  gatherSum(_buffer);

  // for now only run FFT on processor 0
  if (processor_id() != 0)
    return;

  // execute transform
  {
    TIME_SECTION(_perf_fft);
    fftw_execute(_plan);
  }
}

void
FourierTransform::threadJoin(const UserObject & y)
{
  const FourierTransform & vpp = static_cast<const FourierTransform &>(y);
  mooseAssert(vpp._buffer.size() == _buffer.size(), "Inconsistent buffer sizes across threads");

  // sum the two vectors using the binary functional std::plus
  std::transform(vpp._buffer.begin(),
                 vpp._buffer.end(),
                 _buffer.begin(),
                 _buffer.begin(),
                 std::plus<double>());
}

#endif
