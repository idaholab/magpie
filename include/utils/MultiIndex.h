#ifndef MULTIINDEX_H
#define MULTIINDEX_H

#include "Moose.h"
#include "MooseError.h"

namespace MagpieUtils
{

/**
 * Implements a container class for multi-indexed objects
 * with an arbitrary number of indices.
 */
template<class T>
class MultiIndex
{
public:
  MultiIndex(std::vector<unsigned int> shape)
  {
    _shape = shape;
    _dim = _shape.size();

    _nentries = 1;
    for (unsigned int d = 0; d < _dim; ++d)
      _nentries *= _shape[d];

    _data.resize(_nentries);

    // compute the accumulated shapes
    // as[0] = I_1 * I_2 ...* I_{M}, as[1] = I_2 * I_3 ...* I_{M}
    _accumulated_shape.resize(_dim);
    for (unsigned int d = 0; d < _dim; ++d)
    {
      unsigned int k = 1.0;
      for (unsigned int j = d + 1; j < _dim; ++j)
        k *= _shape[j];
      _accumulated_shape[d] = k;
    }
  }

  T & operator() (const std::vector<unsigned int> & indices)
  {
    mooseAssert(indices.size() == _dim, "Indices vector has wrong size. size=" << indices.size() << " vs. dim=" << _dim);
  #if DEBUG
    for (unsigned int j = 0; j < indices.size(); ++j)
      if (indices[j] >= _shape[j])
        mooseError("Indices vector at entry " << j << " is " << indices[j] << " vs. shape " << _shape[j]);
  #endif

    // implement the index
    // index = i_M + i_{M-1} * I_M + i_{M-1} * I_M * I_{M-1} ...
    unsigned int index = 0;
    for (unsigned int d = 0; d < _dim; ++d)
      index += indices[d] * _accumulated_shape[d];

    return _data[index];
  }

  unsigned int size(unsigned int d) const
  {
    mooseAssert(d < _dim, "Argument in size " << d << " larger or equal than " << _dim);
    return _shape[d];
  }

  const std::vector<unsigned int> & shape() const
  {
    return _shape;
  }

  unsigned int dimension() const
  {
    return _dim;
  }

  unsigned int nEntries() const
  {
    return _nentries;
  }

protected:
  /// the number of dimensions
  unsigned int _dim;
  /// the number of entries
  unsigned int _nentries;
  /// the size along each index
  std::vector<unsigned int> _shape;
  /// accumulate index for easier computation of unrolled index
  std::vector<unsigned int> _accumulated_shape;
  /// the data unrolled into a vector
  std::vector<T> _data;
};

} // namespace MagpieUtils

#endif //MULTIINDEX_H
