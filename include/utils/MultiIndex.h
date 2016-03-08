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
  typedef T value_type;
  typedef std::vector<unsigned int> size_type;

  MultiIndex(const size_type & shape);
  MultiIndex(const size_type & shape, const std::vector<T> & data);

  T & operator() (const size_type & indices);
  const T & operator() (const size_type & indices) const;

  // return the container size as dim dimensional vector
  const size_type & size() const { return _shape; }

  unsigned int dim() const { return _shape.size(); }
  unsigned int nEntries() const { return _nentries; }

  // NOTE: resize cannot change the dimensionality of the data.
  // Could possibly implement a reshape function.
  void resize(const size_type & shape);

  // NOTE: assign allows changing the dimensionality of the data.
  void assign(const size_type & shape, T value);

  /**
   * Implement loadHelper and storeHelper for easier data (de)serialization
   */
  void storeMultiIndex(std::ostream & stream, void * context)
  {
    ::dataStore(stream, _shape, context);
    ::dataStore(stream, _data, context);
  }

  void loadMultiIndex(std::istream & stream, void * context)
  {
    ::dataLoad(stream, _shape, context);
    ::dataLoad(stream, _data, context);
    _dim = _shape.size();
    _nentries = _data.size();
    buildAccumulatedShape();
  }

  /**
   * Nested class iterators
   */
  class iterator;

  /**
   * Return iterators for begin and end of this container
   */
  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, _nentries); }

protected:
  /// given a flat index computes the vector of indices i0, i1, i2, ...
  void findIndexVector(unsigned int flat_index, size_type & indices) const
  {
    indices.resize(_dim);
    for (unsigned int d = 0; d < _dim; ++d)
    {
      unsigned int i = flat_index / _accumulated_shape[d];
      indices[d] = i;
      flat_index -= i * _accumulated_shape[d];
    }
  }

  /// the size along each index
  size_type _shape;

  /// the number of dimensions TODO: get rid of this -> _shape.size()
  unsigned int _dim;

  /// the number of entries TODO: get rid of this -> _data.size()
  unsigned int _nentries;

  /// accumulate index for easier computation of unrolled index
  size_type _accumulated_shape;

  /// the data unrolled into a vector
  std::vector<T> _data;

private:
  void buildAccumulatedShape();
  void reshape(const size_type & shape);
};


template <class T>
MultiIndex<T>::MultiIndex(const size_type & shape)
{
  reshape(shape);
  _data.resize(_nentries);
}

template <class T>
MultiIndex<T>::MultiIndex(const size_type & shape, const std::vector<T> & data)
{
  reshape(shape);

  if (data.size() != _nentries)
    mooseError("shape and data arguments' sizes are inconsistent.");
  _data = data;
}

template <class T>
T &
MultiIndex<T>::operator() (const size_type & indices)
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

template <class T>
void
MultiIndex<T>::resize(const size_type & shape)
{
  if (shape.size() != _shape.size())
    mooseError("resize cannot change the dimensionality of MultiIndex.");

  // first copy the old shape and data
  size_type old_shape = _shape;
  size_type old_accumulated_shape = _accumulated_shape;
  std::vector<T> old_data = _data;

  // reset _shape & recompute meta data
  reshape(shape);

  // fill in _data to all possible indices
  _data.assign(_nentries, T(0));
  for (unsigned int j = 0; j < _nentries; ++j)
  {
    size_type indices;
    findIndexVector(j, indices);

    // check if indices existed in the old version of _data
    bool existed_in_old = true;
    for (unsigned int d = 0; d < _dim; ++d)
      if (indices[d] >= old_shape[d])
      {
        existed_in_old = false;
        break;
      }

    // find the corresponding old_j
    if (existed_in_old)
    {
      unsigned int old_j = 0;
      for (unsigned int d = 0; d < _dim; ++d)
        old_j += indices[d] * old_accumulated_shape[d];

      // finally set the data entry
      _data[j] = old_data[old_j];
    }
  }
}

template <class T>
void
MultiIndex<T>::assign(const size_type & shape, T value)
{
  reshape(shape);
  _data.assign(_nentries, value);
}

// compute the accumulated shapes as:
// as[0] = I_1 * I_2 ...* I_{M}, as[1] = I_2 * I_3 ...* I_{M} ...
template <class T>
void
MultiIndex<T>::buildAccumulatedShape()
{
  // TODO: simplify - this is needlessly complicated - can be done in a single loop
  _accumulated_shape.resize(_dim);
  for (unsigned int d = 0; d < _dim; ++d)
  {
    unsigned int k = 1;
    for (unsigned int j = d + 1; j < _dim; ++j)
      k *= _shape[j];
    _accumulated_shape[d] = k;
  }
}

template <class T>
void
MultiIndex<T>::reshape(const size_type & shape)
{
  mooseAssert(shape.size() > 0, "Zero dimensional MultiIndex objects are not supported.");
  _shape = shape;
  _dim = shape.size();

  _nentries = 1;
  for (unsigned int d = 0; d < _dim; ++d)
    _nentries *= _shape[d];

  buildAccumulatedShape();
}

template <class T>
class MultiIndex<T>::iterator
{
public:
  iterator(MultiIndex<T> & multi_index, unsigned int position) :
      _multi_index(multi_index),
      _flat_index(position)
  {
  }

  // Simple data getters
  unsigned int flatIndex() const {return _flat_index;}
  MultiIndex<T> & getMultiIndexObject() const {return _multi_index;}

  // Allow retrieving indices vector from position and single index
  size_type indices() const
  {
    size_type indices;
    _multi_index.findIndexVector(_flat_index, indices);
    return indices;
  }

  unsigned int index(unsigned int d) const
  {
    mooseAssert(d < _multi_index.dimension(), "Dimension d= " << d << " exceeds dim=" << _multi_index.dimension());
    return indices()[d];
  }

  // assignment =
  iterator & operator= (const iterator & other)
  {
    _multi_index = other.getMultiIndexObject();
    _flat_index = other.flatIndex();
    return *this;
  }

  // prefix ++
  iterator & operator++ ()
  {
    ++_flat_index;
    return *this;
  }

  // postfix ++
  iterator & operator++ (int /*k*/)
  {
    iterator clone(*this);
    ++_flat_index;
    return clone;
  }

  // prefix --
  iterator & operator-- ()
  {
    --_flat_index;
    return *this;
  }

  // postfix --
  iterator & operator-- (int /*k*/)
  {
    iterator clone(*this);
    --_flat_index;
    return clone;
  }

  // ==
  // MultiIndexObjects must be exactly the same (pointers to same address)
  // and _flat_index must be the same value
  bool operator== (const iterator & other) const { return _flat_index == other.flatIndex() && &_multi_index == &other.getMultiIndexObject(); }

  // !=
  bool operator!= (const iterator & other) const { return !(*this == other); }

  // overload the dereferencing operator
  T & operator* () { return _multi_index._data[_flat_index]; }

protected:
  MultiIndex<T> & _multi_index;
  unsigned int _flat_index;
};

template<class T>
inline void
dataStore(std::ostream & stream, MultiIndex<T> & ad, void * context) { ad.storeMultiIndex(stream, context); }

template<class T>
inline void
dataLoad(std::istream & stream, MultiIndex<T> & ad, void * context) { ad.loadMultiIndex(stream, context); }

#endif //MULTIINDEX_H
