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

template <>
InputParameters
validParams<FFTBufferAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Fetch values from an FFT buffer.");
  params.addRequiredParam<UserObjectName>("fft_buffer", "Real valued FFT buffer object");
  return params;
}

FFTBufferAux::FFTBufferAux(const InputParameters & parameters)
  : AuxKernel(parameters), _buffer(getUserObject<FFTBufferBase<Real>>("fft_buffer"))
{
}

Real
FFTBufferAux::computeValue()
{
  std::cout << 'B';
  return _buffer(_current_elem->centroid());
}
