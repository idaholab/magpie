/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#include "FourierPowerSpectrum.h"
#include "FourierTransform.h"
#include "MyTRIMMesh.h"

#include "libmesh/utility.h"

#include <functional>

registerMooseObject("MagpieApp", FourierPowerSpectrum);

template <>
InputParameters
validParams<FourierPowerSpectrum>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addClassDescription("Compute the power spectrum of a given fast fourier transform. The "
                             "resulting frequency is in reciprocal length units.");
  params.addRequiredParam<UserObjectName>(
      "fourier_transform", "FourierTransform user object to compute teh power spectrum of");
  return params;
}

FourierPowerSpectrum::FourierPowerSpectrum(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _fourier_transform(getUserObject<FourierTransform>("fourier_transform")),
    _dim(_fourier_transform.dimension()),
    _grid(_fourier_transform.getGrid()),
    _box_size(_fourier_transform.getBoxSize()),
    _buffer(_fourier_transform.getBuffer()),
    _spectrum_size(0),
    _max_frequency(0),
    _power_spectrum(declareVector("power_spectrum")),
    _frequency(declareVector("frequency"))
{
  // estimate spectrum bins
  for (unsigned int i = 0; i < _dim; ++i)
  {
    _max_frequency += Utility::pow<2>(_grid[i] * 0.5 / _box_size(i));
    _spectrum_size += Utility::pow<2>(_grid[i]);
  }

  _max_frequency = std::sqrt(_max_frequency);
  _spectrum_size = std::sqrt(_spectrum_size);

  // resize power spectrum vector
  _power_spectrum.resize(_spectrum_size);
  _frequency.resize(_spectrum_size);
  _n_samples.resize(_spectrum_size);

  // set frequencies
  for (auto i = beginIndex(_frequency); i < _spectrum_size; ++i)
    _frequency[i] = (i * _max_frequency) / _spectrum_size;
}

void
FourierPowerSpectrum::computePowerSpectrum(std::vector<int> & c,
                                           std::vector<Real> & F,
                                           std::size_t i)
{
  if (i == 0)
  {
    Real sum_high = 0.0;
    for (unsigned int j = 1; j < _dim; ++j)
      sum_high += F[j];

    for (c[0] = 0; c[0] < _grid[0]; ++c[0])
    {
      // do stuff at lowest dimension
      const Real sum =
          sum_high + Utility::pow<2>((c[0] * 2 > _grid[0] ? _grid[0] - c[0] : c[0]) / _box_size(0));

      // recalculate index (row major)
      int index = 0;
      for (unsigned int j = 0; j < _dim; ++j)
        index = index * _grid[j] + c[j];

      // find bin and add to spectrum
      const std::size_t bin = std::floor((std::sqrt(sum) * _spectrum_size) / _max_frequency);
      mooseAssert(bin < _power_spectrum.size(), "Bin out of bounds in power spectrum");
      _power_spectrum[bin] += Utility::pow<2>(_buffer[index]);
      ++_n_samples[bin];
    }
  }
  else
    // iterate over lower dimensions
    for (c[i] = 0; c[i] < _grid[i]; ++c[i])
    {
      F[i] = Utility::pow<2>((c[i] * 2 > _grid[i] ? _grid[i] - c[i] : c[i]) / _box_size(i));
      computePowerSpectrum(c, F, i - 1);
    }
}

void
FourierPowerSpectrum::execute()
{
  // clear spectrum
  std::fill(_power_spectrum.begin(), _power_spectrum.end(), 0.0);

  // for now only compute spectrum on processor 0
  if (processor_id() == 0)
  {
    // clear sample count
    std::fill(_n_samples.begin(), _n_samples.end(), 0);

    // compute the power spectrum
    std::vector<int> c(_dim, 0);
    std::vector<Real> F(_dim);
    computePowerSpectrum(c, F, _dim - 1);

    // normalize by sample count
    mooseAssert(_power_spectrum.size() == _n_samples.size(),
                "Size mismatch between spectrum and sample number vectors");
    for (auto i = beginIndex(_power_spectrum); i < _power_spectrum.size(); ++i)
      if (_n_samples[i] > 0)
        _power_spectrum[i] /= _n_samples[i];
  }
}

void
FourierPowerSpectrum::finalize()
{
  // broadcast result to all cores (lazily using gatherSum - all other ranks have zero vectors)
  gatherSum(_power_spectrum);
}

#endif // FFTW3_ENABLED
