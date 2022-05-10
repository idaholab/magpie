/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "FFTBufferAux.h"
#include "FFTBufferBase.h"

registerMooseObject("MagpieApp", FFTBufferAux);

InputParameters
FFTBufferAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Fetch values from an FFT buffer.");
  params.addRequiredParam<UserObjectName>("fft_buffer", "Real valued FFT buffer object");
  params.addParam<std::vector<unsigned int>>(
      "component", std::vector<unsigned int>(), "Component to access for higher order types");
  return params;
}

FFTBufferAux::FFTBufferAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _buffer(getUserObject<UserObject>("fft_buffer")),
    _component(getParam<std::vector<unsigned int>>("component"))
{
  // cast into all possible buffer options
  _real_buffer = dynamic_cast<const FFTBufferBase<Real> *>(&_buffer);
  _realvectorvalue_buffer = dynamic_cast<const FFTBufferBase<RealVectorValue> *>(&_buffer);
  _ranktwotensor_buffer = dynamic_cast<const FFTBufferBase<RankTwoTensor> *>(&_buffer);
  _rankthreetensor_buffer = dynamic_cast<const FFTBufferBase<RankThreeTensor> *>(&_buffer);
  _rankfourtensor_buffer = dynamic_cast<const FFTBufferBase<RankFourTensor> *>(&_buffer);

  // one of those should not result in a nullptr
  if (_real_buffer)
  {
    if (_component.size() != 0)
      paramError("component", "Do not specify a comonent for Real type FFT buffers.");
  }
  else if (_realvectorvalue_buffer)
  {
    if (_component.size() != 1)
      paramError("component",
                 "Specify one index value to access a RealVectorValue type FFT buffer component.");
  }
  else if (_ranktwotensor_buffer)
  {
    if (_component.size() != 2)
      paramError("component",
                 "Specify two index values to access a RankTwoTensor type FFT buffer component.");
  }
  else if (_rankthreetensor_buffer)
  {
    if (_component.size() != 3)
      paramError(
          "component",
          "Specify three index values to access a RankThreeTensor type FFT buffer component.");
  }
  else if (_rankfourtensor_buffer)
  {
    if (_component.size() != 4)
      paramError("component",
                 "Specify four index values to access a RankFourTensor type FFT buffer component.");
  }
  else
    paramError("fft_buffer", "Unsupported buffer type");
}

Real
FFTBufferAux::computeValue()
{
  if (_real_buffer)
    return (*_real_buffer)(_current_elem->vertex_average());
  if (_realvectorvalue_buffer)
    return (*_realvectorvalue_buffer)(_current_elem->vertex_average())(_component[0]);
  if (_ranktwotensor_buffer)
    return (*_ranktwotensor_buffer)(_current_elem->vertex_average())(_component[0], _component[1]);
  if (_rankthreetensor_buffer)
    return (*_rankthreetensor_buffer)(_current_elem->vertex_average())(
        _component[0], _component[1], _component[2]);
  if (_rankfourtensor_buffer)
    return (*_rankfourtensor_buffer)(_current_elem->vertex_average())(
        _component[0], _component[1], _component[2], _component[3]);
  else
    paramError("fft_buffer", "Unsupported buffer type");
}
