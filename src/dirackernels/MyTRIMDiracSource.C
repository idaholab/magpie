#include "MyTRIMDiracSource.h"
#include "MyTRIMDiracRun.h"

template<>
InputParameters validParams<MyTRIMDiracSource>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMDiracRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");

  return params;
}

MyTRIMDiracSource::MyTRIMDiracSource(const InputParameters & parameters) :
    DiracKernel(parameters),
    _mytrim(this->template getUserObject<MyTRIMDiracRun>("runner")),
    _ivar(this->template getParam<unsigned int>("ivar")),
    _defect(this->template getParam<MooseEnum>("defect"))
{
}

void
MyTRIMDiracSource::addPoints()
{
  for (auto && defect : _mytrim.result())
    if (defect._type == _defect && defect._var == _ivar)
      addPoint(defect._location);
}

Real
MyTRIMDiracSource::computeQpResidual()
{
  return -_test[_i][_qp];
}
