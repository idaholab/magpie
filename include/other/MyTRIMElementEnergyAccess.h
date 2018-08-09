/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMELEMENTENERGYACCESS_H
#define MYTRIMELEMENTENERGYACCESS_H

#include "MyTRIMElementRun.h"

/**
 * Interface class ("Veneer") to provide encapsulate fetching energy deposition
 * rates from a MyTRIMElementRun class
 */
template <class T>
class MyTRIMElementEnergyAccess : public T
{
public:
  MyTRIMElementEnergyAccess(const InputParameters & parameters);

  static InputParameters validParams();
  Real getEnergyDensity();

protected:
  const MyTRIMElementRun & _mytrim;

private:
  /// calculate values only for qp 0 and cache them here
  Real _value_cache;
};

template <class T>
MyTRIMElementEnergyAccess<T>::MyTRIMElementEnergyAccess(const InputParameters & parameters)
  : T(parameters), _mytrim(this->template getUserObject<MyTRIMElementRun>("runner"))
{
  if (this->isNodal())
    mooseError("MyTRIMElementEnergyAccess needs to be applied to an elemental AuxVariable.");
}

template <typename T>
InputParameters
MyTRIMElementEnergyAccess<T>::validParams()
{
  InputParameters params = ::validParams<T>();
  params.addRequiredParam<UserObjectName>(
      "runner", "Name of the MyTRIMElementRun userobject to pull data from.");
  return params;
}

template <typename T>
Real
MyTRIMElementEnergyAccess<T>::getEnergyDensity()
{
  if (this->_qp == 0)
  {
    auto & result = _mytrim.result(this->_current_elem);

    // Warning: this volume is only correct in cartesian coordinate systems
    const Real volume = this->_current_elem->volume();

    _value_cache = result._energy / volume;
  }

  return _value_cache;
}

#endif // MYTRIMELEMENTENERGYACCESS_H
