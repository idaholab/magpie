#ifndef MYTRIMRESULTACCESS_H
#define MYTRIMRESULTACCESS_H

#include "MyTRIMRun.h"

/**
 * Interface class ("Veneer") to provide encapsulate fetching defect production
 * rates from a MyTRIMRun class
 */
template <class T>
class MyTRIMResultAccess : public T
{
public:
  MyTRIMResultAccess(const InputParameters & parameters);

  static InputParameters validParams();
  Real getDefectRate();

private:
  const MyTRIMRun & _mytrim;
  const unsigned int _ivar;
  const unsigned int _defect;

  /// calculate values only for qp 0 and cache them here
  Real _value_cache;
};


template <class T>
MyTRIMResultAccess<T>::MyTRIMResultAccess(const InputParameters & parameters) :
    T(parameters),
    _mytrim(getUserObject<MyTRIMRun>("runner")),
    _ivar(getParam<unsigned int>("ivar")),
    _defect(getParam<MooseEnum>("defect"))
{
  if (isNodal())
    mooseError("MyTRIMResultAccess needs to be applied to an elemental AuxVariable.");

  if (_ivar >= _mytrim.nVars())
    mooseError("Requested invalid element index.");
}

template<typename T>
InputParameters
MyTRIMResultAccess<T>::validParams()
{
  InputParameters params = ::validParams<T>();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");
  return params;
}

Real
MyTRIMResultAccess<T>::getDefectRate()
{
  if (_qp == 0)
  {
    const MyTRIMRun::MyTRIMResult & result = _mytrim.result(this->_current_elem);
    mooseAssert(_ivar < result.size(), "Result set does not contain the requested element.");

    const Real volume = this->_current_elem->volume();

    switch (_defect)
    {
      case 0: // vacancy
        _value_cache = result[_ivar].first / volume;
        break;

      case 1: // interstitial
        _value_cache = result[_ivar].second / volume;
        break;

      default:
        mooseError("Internal error.");
    }
  }

  return _value_cache;
}

#endif //MYTRIMRESULTACCESS_H
