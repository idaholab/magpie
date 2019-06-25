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

// Forward Declarations
class IsotopeRecoilRateSampler;
class NeutronicsSpectrumSamplerBase;

template <>
InputParameters validParams<IsotopeRecoilRateSampler>();

class IsotopeRecoilRateSampler : public GeneralVectorPostprocessor
{
public:
  IsotopeRecoilRateSampler(const InputParameters & parameters);

  virtual ~IsotopeRecoilRateSampler() {}

  virtual void initialize() {}
  virtual void execute();

protected:
  std::string _target_isotope;
  std::vector<unsigned int> _point_ids;
  const NeutronicsSpectrumSamplerBase & _neutronics_sampler;
  const PostprocessorValue & _scaling_factor;

  std::vector<Real> & _recoil_rates;
};
