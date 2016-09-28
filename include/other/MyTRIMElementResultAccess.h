#ifndef MYTRIMELEMENTRESULTACCESS_H
#define MYTRIMELEMENTRESULTACCESS_H

#include "MyTRIMElementRun.h"

/**
 * Interface class ("Veneer") to provide encapsulate fetching defect production
 * rates from a MyTRIMElementRun class
 */
template <class T>
class MyTRIMElementResultAccess : public T
{
public:
  MyTRIMElementResultAccess(const InputParameters & parameters);

  static InputParameters validParams();
  Real getDefectRate();

protected:
  const MyTRIMElementRun & _mytrim;
  const unsigned int _ivar;
  const unsigned int _defect;

private:
  /// calculate values only for qp 0 and cache them here
  Real _value_cache;
};


template <class T>
MyTRIMElementResultAccess<T>::MyTRIMElementResultAccess(const InputParameters & parameters) :
    T(parameters),
    _mytrim(this->template getUserObject<MyTRIMElementRun>("runner")),
    _ivar(this->template getParam<unsigned int>("ivar")),
    _defect(this->template getParam<MooseEnum>("defect"))
{
  if (this->isNodal())
    mooseError("MyTRIMElementResultAccess needs to be applied to an elemental AuxVariable.");

  if (_ivar >= _mytrim.nVars())
    mooseError("Requested invalid element index.");
}

template<typename T>
InputParameters
MyTRIMElementResultAccess<T>::validParams()
{
  InputParameters params = ::validParams<T>();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMElementRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");
  return params;
}

template<typename T>
Real
MyTRIMElementResultAccess<T>::getDefectRate()
{
  if (this->_qp == 0)
  {
    auto & result = _mytrim.result(this->_current_elem);
    mooseAssert(_ivar < result._defects.size(), "Result set does not contain the requested element.");

    const Real volume = this->_current_elem->volume();

    switch (_defect)
    {
      case 0: // vacancy
        _value_cache = result._defects[_ivar]._vacancies / volume;
        break;

      case 1: // interstitial
        _value_cache = result._defects[_ivar]._interstitials / volume;
        break;

      default:
        mooseError("Internal error.");
    }
  }

  return _value_cache;
}

#endif //MYTRIMELEMENTRESULTACCESS_H
