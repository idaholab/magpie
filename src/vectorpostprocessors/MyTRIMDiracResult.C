#include "MyTRIMDiracResult.h"
#include "MyTRIMDiracRun.h"

template<>
InputParameters validParams<MyTRIMDiracResult>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMDiracRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");
  return params;
}

MyTRIMDiracResult::MyTRIMDiracResult(const InputParameters & parameters) :
    GeneralVectorPostprocessor(parameters),
    _mytrim(getUserObject<MyTRIMDiracRun>("runner")),
    _ivar(getParam<unsigned int>("ivar")),
    _defect(getParam<MooseEnum>("defect")),
    _x(declareVector("x")),
    _y(declareVector("y")),
    _z(declareVector("z"))
{
}

void
MyTRIMDiracResult::initialize()
{
  _x.clear();
  _y.clear();
  _z.clear();
}

void
MyTRIMDiracResult::execute()
{
  for (auto && defect : _mytrim.result())
    if (defect._type == _defect && defect._var == _ivar)
    {
      _x.push_back(defect._location(0));
      _y.push_back(defect._location(1));
      _z.push_back(defect._location(2));
    }
}

void
MyTRIMDiracResult::finalize()
{
  // the MyTRIMDiracRunner already does the necessary communication
}

void
MyTRIMDiracResult::threadJoin(const UserObject &)
{
  // the MyTRIMDiracRunner already does the necessary communication
}
