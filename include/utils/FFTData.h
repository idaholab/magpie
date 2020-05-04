/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ComplexTypes.h"

// helper template to select the corresponding scalar type
template <typename T>
struct FFTScalarType
{
  using type = Real;
};

template <typename T>
struct FFTScalarType<std::complex<T>>
{
  using type = Complex;
};

/**
 * Helper class to hold the reciprocal space data
 */
template <typename T>
class FFTData
{
  using ScalarT = typename FFTScalarType<T>::type;

public:
  FFTData(std::size_t size = 0) { resize(size); }
  void resize(std::size_t size) { _buffer.resize(size); }

  ///@{ data access by index
  const T & operator[](std::size_t i) const { return _buffer[i]; }
  T & operator[](std::size_t i) { return _buffer[i]; }
  ///@}

  ///@{ convenience math operators
  FFTData<T> & operator+=(FFTData<T> const & rhs);
  FFTData<T> & operator-=(FFTData<T> const & rhs);
  FFTData<T> & operator*=(FFTData<ScalarT> const & rhs);
  FFTData<T> & operator/=(FFTData<ScalarT> const & rhs);
  FFTData<T> & operator*=(Real rhs);
  FFTData<T> & operator/=(Real rhs);
  FFTData<T> & operator=(FFTData<T> const & rhs);
  ///@}

  /// return the number of proper grid cells
  std::size_t size() const { return _buffer.size(); }

  /// get the addres of the first data element of the ith object in the buffer
  void * start(std::size_t i);

  /// get the number of transforms required for type T
  std::size_t howMany() const;

protected:
  /// FFT data buffer
  std::vector<T> _buffer;
};
