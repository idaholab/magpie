/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#ifndef FOURIERPOWERSPECTRUM_H
#define FOURIERPOWERSPECTRUM_H

#include "GeneralVectorPostprocessor.h"

#include <memory>

class FourierPowerSpectrum;
class FourierTransform;

template <>
InputParameters validParams<FourierPowerSpectrum>();

/**
 * Compute the power spectrum from the data of a FourierTransform object.
 */
class FourierPowerSpectrum : public GeneralVectorPostprocessor
{
public:
  FourierPowerSpectrum(const InputParameters & parameters);

  virtual void initialize() override{};
  virtual void execute() override;
  virtual void finalize() override;

protected:
  void computePowerSpectrum(std::vector<int> & c, std::vector<Real> & F, std::size_t i);

  /// Fourier Transform to provide the data
  const FourierTransform & _fourier_transform;

  /// mesh dimension
  const unsigned int _dim;

  /// grid size for FFT (needs to be signed for FFTW)
  const std::vector<int> & _grid;

  /// simulation box extents
  const Point _box_size;

  /// FFTW data buffer
  const std::vector<double> & _buffer;

  ///{@ power spectrum data
  std::size_t _spectrum_size;
  Real _max_frequency;
  VectorPostprocessorValue & _power_spectrum;
  VectorPostprocessorValue & _frequency;
  std::vector<int> _n_samples;
  ///@}
};

#endif // FOURIERPOWERSPECTRUM_H
#endif // FFTW3_ENABLED
