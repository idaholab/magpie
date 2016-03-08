#ifndef MULTIINDEX_H
#define MULTIINDEX_H

#include "Moose.h"
#include "MooseError.h"
#include "DataIO.h"

/**
 * Implements a container class for multi-indexed objects
 * with an arbitrary number of indices.
 */
template<class T>
class MultiIndex
{
public:
  typedef unsigned int size_type;

  MultiIndex(std::vector<size_type> shape) :
    _shape(shape),
    _dim(_shape.size())
  {
    _nentries = 1;
    for (size_type d = 0; d < _dim; ++d)
      _nentries *= _shape[d];

    _data.resize(_nentries);

    // compute the accumulated shapes as:
    // as[0] = I_1 * I_2 ...* I_{M}, as[1] = I_2 * I_3 ...* I_{M} ...
    _accumulated_shape.resize(_dim);
    for (size_type d = 0; d < _dim; ++d)
    {
      size_type k = 1.0;
      for (size_type j = d + 1; j < _dim; ++j)
        k *= _shape[j];
      _accumulated_shape[d] = k;
    }
  }

  MultiIndex(std::vector<size_type> shape, std::vector<T> data) :
    _shape(shape),
    _dim(_shape.size())
  {
    _nentries = 1;
    for (size_type d = 0; d < _dim; ++d)
      _nentries *= _shape[d];

    if (data.size() != _nentries)
      mooseError("shape and data arguments' sizes are inconsistent.");
    _data = data;

    // compute the accumulated shapes as:
    // as[0] = I_1 * I_2 ...* I_{M}, as[1] = I_2 * I_3 ...* I_{M} ...
    _accumulated_shape.resize(_dim);
    for (size_type d = 0; d < _dim; ++d)
    {
      size_type k = 1;
      for (size_type j = d + 1; j < _dim; ++j)
        k *= _shape[j];
      _accumulated_shape[d] = k;
    }
  }

  T & operator() (const std::vector<size_type> & indices)
  {
    mooseAssert(indices.size() == _dim, "Indices vector has wrong size. size=" << indices.size() << " vs. dim=" << _dim);
#if DEBUG
    for (size_type j = 0; j < indices.size(); ++j)
      if (indices[j] >= _shape[j])
        mooseError("Indices vector at entry " << j << " is " << indices[j] << " vs. shape " << _shape[j]);
#endif

    // implement the index
    // index = i_M + i_{M-1} * I_M + i_{M-1} * I_M * I_{M-1} ...
    size_type index = 0;
    for (size_type d = 0; d < _dim; ++d)
      index += indices[d] * _accumulated_shape[d];

    return _data[index];
  }

  size_type size(size_type d) const
  {
    mooseAssert(d < _dim, "Argument in size " << d << " larger or equal than " << _dim);
    return _shape[d];
  }

  const std::vector<size_type> & shape() const { return _shape; }

  size_type dimension() const { return _dim; }

  size_type nEntries() const { return _nentries; }

  // NOTE: resize cannot change the dimensionality of the data.
  // Could possibly implement a reshape function.
  void resize(std::vector<size_type> shape)
  {
    if (shape.size() != _dim)
      mooseError("resize cannot change the dimensionality of MultiIndex.");

    // first copy the old shape and data
    std::vector<size_type> old_shape = _shape;
    std::vector<size_type> old_accumulated_shape = _accumulated_shape;
    std::vector<T> old_data = _data;

    // reset _shape & recompute meta data
    _shape = shape;
    _nentries = 1;
    for (size_type d = 0; d < _dim; ++d)
      _nentries *= _shape[d];

    for (size_type d = 0; d < _dim; ++d)
    {
      size_type k = 1.0;
      for (size_type j = d + 1; j < _dim; ++j)
        k *= _shape[j];
      _accumulated_shape[d] = k;
    }

    // fill in _data to all possible indices
    _data.assign(_nentries, 0.0);
    for (size_type j = 0; j < _nentries; ++j)
    {
      std::vector<size_type> indices;
      findIndexVector(j, indices);

      // check if indices existed in the old version of _data
      bool existed_in_old = true;
      for (size_type d = 0; d < _dim; ++d)
        if (indices[d] >= old_shape[d])
        {
          existed_in_old = false;
          break;
        }

      // find the corresponding old_j
      if (existed_in_old)
      {
        size_type old_j = 0;
        for (size_type d = 0; d < _dim; ++d)
          old_j += indices[d] * old_accumulated_shape[d];

        // finally set the data entry
        _data[j] = old_data[old_j];
      }
    }
  }

  // NOTE: assign allows changing the dimensionality of the data.
  void assign(std::vector<size_type> shape, T value)
  {
    _shape = shape;
    _dim = _shape.size();

    _nentries = 1;
    for (size_type d = 0; d < _dim; ++d)
      _nentries *= _shape[d];

    _data.assign(_nentries, value);

    // compute the accumulated shapes as:
    // as[0] = I_1 * I_2 ...* I_{M}, as[1] = I_2 * I_3 ...* I_{M} ...
    _accumulated_shape.resize(_dim);
    for (size_type d = 0; d < _dim; ++d)
    {
      size_type k = 1.0;
      for (size_type j = d + 1; j < _dim; ++j)
        k *= _shape[j];
      _accumulated_shape[d] = k;
    }
  }

  /**
   * Implement loadHelper and storeHelper for easier data (de)serialization
   */
  void storeMultiIndex(std::ostream & stream, void * context)
  {
    storeHelper<size_type>(stream, _dim, context);
    storeHelper<size_type>(stream, _nentries, context);
    storeHelper<std::vector<size_type> >(stream, _shape, context);
    storeHelper<std::vector<size_type> >(stream, _accumulated_shape, context);
    storeHelper<std::vector<T> >(stream, _data, context);
  }

  void loadMultiIndex(std::istream & stream, void * context)
  {
    loadHelper<size_type>(stream, _dim, context);
    loadHelper<size_type>(stream, _nentries, context);
    loadHelper<std::vector<size_type> >(stream, _shape, context);
    loadHelper<std::vector<size_type> >(stream, _accumulated_shape, context);
    loadHelper<std::vector<T> >(stream, _data, context);
  }

  /**
   * Nested class iterators
   */
  class iterator
  {
  public:
    iterator(MultiIndex<T> & multi_index, size_type position) :
     _multi_index(multi_index),
     _current_pos(position)
    {
    }

    // Simple data getters
    size_type getCurrentPos() const {return _current_pos;}
    MultiIndex<T> & getMultiIndexObject() const {return _multi_index;}

    // Allow retrieving indices vector from position and single index
    std::vector<size_type> indices() const
    {
      std::vector<size_type> indices;
      _multi_index.findIndexVector(_current_pos, indices);
      return indices;
    }

    size_type index(size_type d) const
    {
      mooseAssert(d < _multi_index.dimension(), "Dimension d= " << d << " exceeds dim=" << _multi_index.dimension());
      std::vector<size_type> indices;
      _multi_index.findIndexVector(_current_pos, indices);
      return indices[d];
    }

    /*
     * overload the =, ++, --, ==, != operators
     */

    // assignment =
    iterator & operator= (const iterator & other)
    {
      _multi_index = other.getMultiIndexObject();
      _current_pos = other.getCurrentPos();
      return *this;
    }

    // prefix ++
    iterator & operator++ ()
    {
      ++_current_pos;
      return *this;
    }

    // postfix ++
    iterator & operator++ (int /*k*/)
    {
      iterator clone(*this);
      ++_current_pos;
      return clone;
    }

    // prefix --
    iterator & operator-- ()
    {
      --_current_pos;
      return *this;
    }

    // postfix --
    iterator & operator-- (int /*k*/)
    {
      iterator clone(*this);
      --_current_pos;
      return clone;
    }

    // ==
    // MultiIndexObjects must be exactly the same (pointers to same address)
    // and _current_pos must be the same value
    bool operator== (const iterator & other) const { return _current_pos == other.getCurrentPos() && &_multi_index == &other.getMultiIndexObject(); }

    // !=
    bool operator!= (const iterator & other) const { return !(*this == other); }

    // overload the dereferencing operator
    T & operator* () { return _multi_index._data[_current_pos]; }

  protected:
    MultiIndex<T> & _multi_index;
    size_type _current_pos;
   };

  /**
   * Return iterators for begin and end of this container
   */
  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, _nentries); }

protected:
  /// given a flat index computes the vector of indices i0, i1, i2, ...
  void findIndexVector(size_type flat_index, std::vector<size_type> & indices) const
  {
    indices.resize(_dim);
    for (size_type d = 0; d < _dim; ++d)
    {
      size_type i = flat_index / _accumulated_shape[d];
      indices[d] = i;
      flat_index -= i * _accumulated_shape[d];
    }
  }

  /// the size along each index
  std::vector<size_type> _shape;

  /// the number of dimensions
  size_type _dim;

  /// the number of entries
  size_type _nentries;

  /// accumulate index for easier computation of unrolled index
  std::vector<size_type> _accumulated_shape;

  /// the data unrolled into a vector
  std::vector<T> _data;
};

template<class T>
inline void
dataStore(std::ostream & stream, MultiIndex<T> & ad, void * context) { ad.storeMultiIndex(stream, context); }

template<class T>
inline void
dataLoad(std::istream & stream, MultiIndex<T> & ad, void * context) { ad.loadMultiIndex(stream, context); }

#endif //MULTIINDEX_H
