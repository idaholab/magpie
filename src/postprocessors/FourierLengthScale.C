/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#include "FourierLengthScale.h"
#include "FourierTransform.h"
#include "MyTRIMMesh.h"

#include "libmesh/utility.h"

#include <functional>

registerMooseObject("MagpieApp", FourierLengthScale);

InputParameters
FourierLengthScale::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Compute the average length scale form a given Fourier transform in length units.");
  params.addRequiredParam<UserObjectName>(
      "fourier_transform", "FourierTransform user object to compute the length scale of");
  return params;
}

FourierLengthScale::FourierLengthScale(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _fourier_transform(getUserObject<FourierTransform>("fourier_transform")),
    _dim(_fourier_transform.dimension()),
    _grid(_fourier_transform.getGrid()),
    _box_size(_fourier_transform.getBoxSize()),
    _buffer(_fourier_transform.getBuffer())
{
}

void
FourierLengthScale::computeLengthScale(std::vector<int> & c, std::vector<Real> & F, std::size_t i)
{
  if (i == _dim - 1)
  {
    Real sum_high = 0.0;
    for (unsigned int j = 1; j < _dim; ++j)
      sum_high += F[j];

    auto & cc = c[_dim - 1];
    for (cc = 0; cc < _grid[_dim - 1]; ++cc)
    {
      // do stuff at highest dimension
      F[_dim - 1] =
          Utility::pow<2>((cc * 2 > _grid[_dim - 1] ? _grid[_dim - 1] - cc : cc) / _box_size(0));
      Real sum = sum_high + F[_dim - 1];

      // find bin and add to spectrum
      if (sum > 0)
      {
        // compute frequency
        const Real frequency = std::sqrt(sum);

        // compute number of wavevectors contributing to this frequency (frequency^(_dim-1))
        Real normalization = 1.0;
        for (unsigned int j = 1; j < _dim; ++j)
          normalization *= frequency;

        // weight is the normalized intensity at this frequency
        const Real weight = Utility::pow<2>(_buffer[_index]) / normalization;

        // compute weighted average of frequencies
        _length_scale += frequency * weight;
        _weight_sum += weight;
      }
      _index++;
    }
  }
  else
    // iterate over lower dimensions
    for (c[i] = 0; c[i] < _grid[i]; ++c[i])
    {
      F[i] = Utility::pow<2>((c[i] * 2 > _grid[i] ? _grid[i] - c[i] : c[i]) / _box_size(i));
      computeLengthScale(c, F, i + 1);
    }
}

void
FourierLengthScale::execute()
{
  // clear data
  _length_scale = 0.0;
  _weight_sum = 0.0;

  // for now only compute length scale on processor 0
  if (processor_id() == 0)
  {
    // compute the length scale
    std::vector<int> c(_dim, 0);
    std::vector<Real> F(_dim);
    _index = 0;
    computeLengthScale(c, F, 0);

    // divide by weight and take reciprocal (frequency -> length)
    _length_scale = _weight_sum / _length_scale;
  }
}

void
FourierLengthScale::finalize()
{
  // broadcast result to all cores (lazily using gatherSum - all other ranks have zero length scale)
  gatherSum(_length_scale);
}

#endif // FFTW3_ENABLED
