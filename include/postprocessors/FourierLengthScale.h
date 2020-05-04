/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#pragma once

#include "GeneralPostprocessor.h"

#include <memory>

class FourierTransform;

/**
 * Compute the average length scale form a fast fourier transform
 */
class FourierLengthScale : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  FourierLengthScale(const InputParameters & parameters);

  virtual void initialize() override{};
  virtual void execute() override;
  virtual void finalize() override;

  virtual PostprocessorValue getValue() override { return _length_scale; }

protected:
  void computeLengthScale(std::vector<int> & c, std::vector<Real> & F, std::size_t i);

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

  /// length scale computed from the fourier transform
  PostprocessorValue _length_scale;

  /// averaging weight
  Real _weight_sum;

  /// index variable for recursive buffer traversal
  std::size_t _index;
};

#endif // FFTW3_ENABLED
