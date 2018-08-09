/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef ISOTOPERECOILRATE_H
#define ISOTOPERECOILRATE_H

#include "GeneralPostprocessor.h"

// forward declarations
class IsotopeRecoilRate;
class NeutronicsSpectrumSamplerBase;

template <>
InputParameters validParams<IsotopeRecoilRate>();

class IsotopeRecoilRate : public GeneralPostprocessor
{
public:
  IsotopeRecoilRate(const InputParameters & parameters);
  virtual void execute() override {}
  virtual void initialize() override {}
  virtual Real getValue() override;

protected:
  std::string _target_isotope;
  unsigned int _point_id;
  const NeutronicsSpectrumSamplerBase & _neutronics_sampler;
  const PostprocessorValue & _scaling_factor;
};

#endif // IsotopeRecoilRate_H
