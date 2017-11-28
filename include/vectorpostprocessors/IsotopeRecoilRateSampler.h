/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef ISOTOPERECOILRATESAMPLER_H
#define ISOTOPERECOILRATESAMPLER_H

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

  std::vector<Real> _recoil_rates;
};

#endif
