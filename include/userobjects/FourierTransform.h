/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef FFTW3_ENABLED

#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include "ElementUserObject.h"
#include "fftw3.h"

#include <memory>

class FourierTransform;

template <>
InputParameters validParams<FourierTransform>();

/**
 * Compute the fourier transform of a selected variable field.
 */
class FourierTransform : public ElementUserObject
{
public:
  FourierTransform(const InputParameters & parameters);
  ~FourierTransform();

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

  ///@{ public API
  unsigned int dimension() const { return _dim; }
  const std::vector<int> & getGrid() const { return _grid; }
  const Point & getBoxSize() const { return _box_size; }
  const std::vector<double> & getBuffer() const { return _buffer; }
  ///@}

protected:
  /// variable to compute FFT of
  const VariableValue & _var;

  /// mesh dimension
  unsigned int _dim;

  /// grid size for FFT (needs to be signed for FFTW)
  std::vector<int> _grid;

  ///@{ simulation box extents
  Point _min_corner;
  Point _max_corner;
  Point _box_size;
  ///@}

  /// FFT grid cell volume
  Real _cell_volume;

  /// FFTW data buffer
  std::vector<double> _buffer;
  std::size_t _buffer_size;

  /// FFTW plan
  fftw_plan _plan;

  ///@{ timers
  PerfID _perf_plan;
  PerfID _perf_fft;
  ///@}
};

#endif // FOURIERTRANSFORM_H
#endif // FFTW3_ENABLED
