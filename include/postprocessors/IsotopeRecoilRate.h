/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralPostprocessor.h"

class NeutronicsSpectrumSamplerBase;

class IsotopeRecoilRate : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  IsotopeRecoilRate(const InputParameters & parameters);

  virtual void execute() override {}
  virtual void initialize() override {}

  using Postprocessor::getValue;
  virtual Real getValue() const override;

protected:
  std::string _target_isotope;
  unsigned int _point_id;
  const NeutronicsSpectrumSamplerBase & _neutronics_sampler;
  const PostprocessorValue & _scaling_factor;
};
