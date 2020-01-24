/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ElementUserObject.h"

template <typename T>
class FFTBufferBase;

#define usingFFTBufferBaseMembers                                                                  \
  using ElementUserObject::_perf_graph;                                                            \
  using FFTBufferBase<T>::_dim;                                                                    \
  using FFTBufferBase<T>::_grid;                                                                   \
  using FFTBufferBase<T>::_buffer;                                                                 \
  using FFTBufferBase<T>::_start;                                                                  \
  using FFTBufferBase<T>::_stride;                                                                 \
  using FFTBufferBase<T>::_how_many

/**
 * Generic FFT interleaved data buffer base class
 */
template <typename T>
class FFTBufferBase : public ElementUserObject
{
public:
  static InputParameters validParams();

  FFTBufferBase(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute();
  virtual void finalize() {}
  virtual void threadJoin(const UserObject &) {}

  ///@{ transforms
  virtual void forward() = 0;
  virtual void backward() = 0;
  ///@}

  ///@{ data access by index
  const T & operator[](std::size_t i) const { return _buffer[i]; }
  T & operator[](std::size_t i) { return _buffer[i]; }
  ///@}

  ///@{ data access by location
  const T & operator()(const Point & p) const;
  T & operator()(const Point & p);
  ///@}

protected:
  /// get the addres of the first data element of the ith object in the bufefr
  Real * start(std::size_t i);

  /// get the number of transforms required for type T
  std::size_t howMany() const;

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

  ///@{ FFT data buffer
  std::vector<T> _buffer;
  std::size_t _buffer_size;
  ///@}

  /// pointer to the start of the data
  Real * _start;

  /// stride in units of double size
  std::ptrdiff_t _stride;

  /// optional moose sister variabe (to obtain IC from)
  std::vector<const VariableValue *> _moose_variable;

  /// cache the howMany value
  const std::size_t _how_many;
};
