/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralUserObject.h"

template <typename T>
class FFTBufferBase;

/**
 * Generic FFT interleaved data buffer base class
 */
template <typename T>
class FFTBufferBase : public GeneralUserObject
{
public:
  static InputParameters validParams();

  FFTBufferBase(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}

  // transforms
  virtual void forward() = 0;
  virtual void backward() = 0;

  // data access
  const T & operator[](std::size_t i) const { return _buffer[i]; }
  T & operator[](std::size_t i) { return _buffer[i]; }

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
};
