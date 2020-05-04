/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "FFTData.h"

template <typename T>
FFTData<T> &
FFTData<T>::operator+=(FFTData<T> const & rhs)
{
  for (std::size_t i = 0; i < _buffer.size(); ++i)
    _buffer[i] += rhs[i];
  return *this;
}

template <typename T>
FFTData<T> &
FFTData<T>::operator-=(FFTData<T> const & rhs)
{
  for (std::size_t i = 0; i < _buffer.size(); ++i)
    _buffer[i] -= rhs[i];
  return *this;
}

template <typename T>
FFTData<T> &
FFTData<T>::operator*=(FFTData<typename FFTScalarType<T>::type> const & rhs)
{
  for (std::size_t i = 0; i < _buffer.size(); ++i)
    _buffer[i] *= rhs[i];
  return *this;
}

template <typename T>
FFTData<T> &
FFTData<T>::operator/=(FFTData<typename FFTScalarType<T>::type> const & rhs)
{
  for (std::size_t i = 0; i < _buffer.size(); ++i)
    _buffer[i] /= rhs[i];
  return *this;
}

template <typename T>
FFTData<T> &
FFTData<T>::operator*=(Real rhs)
{
  for (std::size_t i = 0; i < _buffer.size(); ++i)
    _buffer[i] *= rhs;
  return *this;
}

template <typename T>
FFTData<T> &
FFTData<T>::operator/=(Real rhs)
{
  const Real reciprocal = 1 / rhs;
  for (std::size_t i = 0; i < _buffer.size(); ++i)
    _buffer[i] *= reciprocal;
  return *this;
}

template <typename T>
FFTData<T> &
FFTData<T>::operator=(FFTData<T> const & rhs)
{
  mooseAssert(size() == rhs.size(), "Buffers need to have identical size!");
  _buffer = rhs._buffer;
  return *this;
}

template <>
void *
FFTData<Real>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i]);
}
template <>
void *
FFTData<Complex>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i]);
}

template <>
void *
FFTData<RealVectorValue>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0));
}
template <>
void *
FFTData<ComplexVectorValue>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0));
}

template <>
void *
FFTData<RankTwoTensor>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0, 0));
}
template <>
void *
FFTData<ComplexRankTwoTensor>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0, 0));
}

template <>
void *
FFTData<RankThreeTensor>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0, 0, 0));
}
template <>
void *
FFTData<ComplexRankThreeTensor>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0, 0, 0));
}

template <>
void *
FFTData<RankFourTensor>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0, 0, 0, 0));
}
template <>
void *
FFTData<ComplexRankFourTensor>::start(std::size_t i)
{
  return reinterpret_cast<void *>(&_buffer[i](0, 0, 0, 0));
}

template <typename T>
std::size_t
FFTData<T>::howMany() const
{
  mooseError("Only call howMany() on the real space buffer!");
}

template <>
std::size_t
FFTData<Real>::howMany() const
{
  return 1;
}

template <>
std::size_t
FFTData<RealVectorValue>::howMany() const
{
  return LIBMESH_DIM;
}

template <>
std::size_t
FFTData<RankTwoTensor>::howMany() const
{
  return LIBMESH_DIM * LIBMESH_DIM;
}

template <>
std::size_t
FFTData<RankThreeTensor>::howMany() const
{
  return LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM;
}

template <>
std::size_t
FFTData<RankFourTensor>::howMany() const
{
  return LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM * LIBMESH_DIM;
}

// explicit instantiations
template class FFTData<Real>;
template class FFTData<RealVectorValue>;
template class FFTData<RankTwoTensor>;
template class FFTData<RankThreeTensor>;
template class FFTData<RankFourTensor>;

template class FFTData<Complex>;
template class FFTData<ComplexVectorValue>;
template class FFTData<ComplexRankTwoTensor>;
template class FFTData<ComplexRankThreeTensor>;
template class FFTData<ComplexRankFourTensor>;
