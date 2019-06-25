/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMELEMENTRESULTACCESS_H
#define MYTRIMELEMENTRESULTACCESS_H

#include "MyTRIMElementRun.h"
#include "MyTRIMRasterizer.h"
#include "ThreadedRecoilLoopBase.h"

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
  const MyTRIMRasterizer & _rasterizer;
  const unsigned int _ivar;
  ThreadedRecoilLoopBase::DefectType _defect;

private:
  /// calculate values only for qp 0 and cache them here
  Real _value_cache;
};

template <class T>
MyTRIMElementResultAccess<T>::MyTRIMElementResultAccess(const InputParameters & parameters)
  : T(parameters),
    _mytrim(getUserObject<MyTRIMElementRun>("runner")),
    _rasterizer(_mytrim.rasterizer()),
    _ivar(getParam<unsigned int>("ivar")),
    _defect(getParam<MooseEnum>("defect").template getEnum<ThreadedRecoilLoopBase::DefectType>())
{
  if (this->isNodal())
    mooseError("MyTRIMElementResultAccess needs to be applied to an elemental AuxVariable.");

  if (_ivar >= _mytrim.nVars())
    mooseError("Requested invalid element index.");
}

template <typename T>
InputParameters
MyTRIMElementResultAccess<T>::validParams()
{
  InputParameters params = ::validParams<T>();
  params.addRequiredParam<UserObjectName>(
      "runner", "Name of the MyTRIMElementRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC=0 INT REPLACEMENT_IN REPLACEMENT_OUT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");
  return params;
}

template <typename T>
Real
MyTRIMElementResultAccess<T>::getDefectRate()
{
  if (this->_qp == 0)
  {
    auto & result = _mytrim.result(this->_current_elem);
    mooseAssert(_ivar < result._defects.size(),
                "Result set does not contain the requested element.");

    const Real volume = this->_current_elem->volume();
    _value_cache = result._defects[_ivar][_defect] / volume;
  }

  return _value_cache;
}

#endif // MYTRIMELEMENTRESULTACCESS_H
