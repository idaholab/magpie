/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

// MOOSE includes
#include "GeneralVectorPostprocessor.h"

class NeutronicsSpectrumSamplerBase;

class IsotopeRecoilRateSampler : public GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  IsotopeRecoilRateSampler(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute();

protected:
  std::string _target_isotope;
  std::vector<unsigned int> _point_ids;
  const NeutronicsSpectrumSamplerBase & _neutronics_sampler;
  const PostprocessorValue & _scaling_factor;

  std::vector<Real> & _recoil_rates;
};
