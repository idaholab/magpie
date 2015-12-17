#include "MyTRIMAux.h"
#include "MyTRIMRun.h"

template<>
InputParameters validParams<MyTRIMAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("runner", "Name of the MyTRIMRun userobject to pull data from.");
  params.addParam<unsigned int>("ivar", "Element index");
  MooseEnum defectType("VAC INT", "VAC");
  params.addParam<MooseEnum>("defect", defectType, "Defect type to read out");
  return params;
}

MyTRIMAux::MyTRIMAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _mytrim(getUserObject<MyTRIMRun>("runner")),
    _ivar(getParam<unsigned int>("ivar")),
    _defect(getParam<MooseEnum>("defect"))
{
  if (isNodal())
    mooseError("MyTRIMAux needs to be applied to an elemental AuxVariable.");

  if (_ivar >= _mytrim.nVars())
    mooseError("Requested invalid element index.");
}

Real
MyTRIMAux::computeValue()
{
  if (_qp == 0)
  {
    const MyTRIMRun::MyTRIMResult & result = _mytrim.result(_current_elem);
    mooseAssert(_ivar < result.size(), "Result set does not contain the requested element.");

    switch(_defect)
    {
      case 0: // vacancy
        _value_cache = result[_ivar].first;
        break;

      case 1: // interstitial
        _value_cache = result[_ivar].second;
        break;

      default:
        mooseError("Internal error.");
    }
  }

  return _value_cache;
}
