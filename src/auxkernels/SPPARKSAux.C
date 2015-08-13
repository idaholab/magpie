#include "SPPARKSAux.h"
#include "SPPARKSUserObject.h"

template<>
InputParameters validParams<SPPARKSAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("user_object", "Name of SPPARKSUserObject");
  params.addParam<unsigned int>("ivar", "Index into SPPARKS iarray");
  params.addParam<unsigned int>("dvar", "Index into SPPARKS darray");
  return params;
}

SPPARKSAux::SPPARKSAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _spparks(getUserObject<SPPARKSUserObject>("user_object")),
    _ivar(isParamValid("ivar") ? getParam<unsigned>("ivar") : -1),
    _dvar(isParamValid("dvar") ? getParam<unsigned>("dvar") : -1)
{
  if ((_ivar >= 0 && _dvar >=0) || (_ivar < 0 && _dvar < 0))
    mooseError("Error in SPPARKSAux, " << _name << ": Either ivar or dvar must be given.");
}

Real
SPPARKSAux::computeValue()
{
  Real value = 0;

  if (_ivar > -1)
    value = _spparks.getIntValue(_current_node->id(), _ivar);
  else
    value = _spparks.getDoubleValue(_current_node->id(), _dvar);

  return value;
}
