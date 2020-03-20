/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ComplexTypes.h"
#include "ElementUserObject.h"
#include "FFTData.h"

template <typename T>
class FFTBufferBase;

#define usingFFTBufferBaseMembers                                                                  \
  using ElementUserObject::_perf_graph;                                                            \
  using FFTBufferBase<T>::_dim;                                                                    \
  using FFTBufferBase<T>::_grid;                                                                   \
  using FFTBufferBase<T>::_real_space_data;                                                        \
  using FFTBufferBase<T>::_reciprocal_space_data;                                                  \
  using FFTBufferBase<T>::_real_space_data_start;                                                  \
  using FFTBufferBase<T>::_reciprocal_space_data_start;                                            \
  using FFTBufferBase<T>::_real_space_data_stride;                                                 \
  using FFTBufferBase<T>::_reciprocal_space_data_stride;                                           \
  using FFTBufferBase<T>::_how_many

/**
 * Generic FFT interleaved data buffer base class
 */
template <typename T>
class FFTBufferBase : public ElementUserObject
{
  using ComplexT = typename ComplexType<T>::type;

public:
  static InputParameters validParams();

  FFTBufferBase(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute();
  virtual void finalize() {}
  virtual void threadJoin(const UserObject &) {}

  ///@{ transforms with proper scaling guaranteed
  virtual void forward();
  virtual void backward();
  ///@}

  ///@{ transforms without proper scaling guaranteed
  virtual void forwardRaw() = 0;
  virtual void backwardRaw() = 0;
  ///@}

  ///@{ scaling required after respective transform
  virtual Real forwardScale() { return 1.0; }
  virtual Real backwardScale() { return forwardScale(); }
  ///@}

  ///@{ buffer access
  FFTData<T> & realSpace() { return _real_space_data; }
  FFTData<ComplexT> & reciprocalSpace() { return _reciprocal_space_data; }
  const FFTData<T> & realSpace() const { return _real_space_data; }
  const FFTData<ComplexT> & reciprocalSpace() const { return _reciprocal_space_data; }

  ///@{ real space data access by location
  const T & operator()(const Point & p) const;
  T & operator()(const Point & p);
  ///@}

  /// return the number of grid cells along each dimension without padding
  const std::vector<int> & grid() const { return _grid; }

  /// return the buffer dimension
  const unsigned int & dim() const { return _dim; }

  /// return the buffer dimension
  const std::vector<Real> & kTable(unsigned int d) const
  {
    mooseAssert(d < _dim, "invalid kTable component");
    return _k_table[d];
  }

protected:
  ///@{ mesh data
  MooseMesh & _mesh;
  unsigned int _dim;
  ///@}

  /// grid size for FFT (needs to be signed for FFTW)
  std::vector<int> _grid;

  ///@{ simulation box extents
  Point _min_corner;
  Point _max_corner;
  Point _box_size;
  ///@}

  /// FFT grid cell volume
  Real _cell_volume;

  ///@{ FFT data buffer and unpadded number of grid cells
  FFTData<T> _real_space_data;
  FFTData<ComplexT> _reciprocal_space_data;
  ///@}

  /// pointer to the start of the data
  Real * _real_space_data_start;
  Complex * _reciprocal_space_data_start;

  /// stride in units of double size
  std::ptrdiff_t _real_space_data_stride;
  std::ptrdiff_t _reciprocal_space_data_stride;

  /// optional moose sister variabe (to obtain IC from)
  std::vector<const VariableValue *> _moose_variable;

  /// pretabulated k-vector components
  std::vector<std::vector<Real>> _k_table;

  /// cache the howMany value
  const std::size_t _how_many;
};
