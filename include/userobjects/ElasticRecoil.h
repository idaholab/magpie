/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/
#ifdef GSL_ENABLED

#pragma once

#include "ScatteringRecoilCrossSection.h"

class ElasticRecoil : public ScatteringRecoilCrossSection
{
public:
  static InputParameters validParams();

  ElasticRecoil(const InputParameters & parameters);

  virtual void execute() override;

protected:
  /// helper function to get mu_C from T and E
  virtual Real getCMCosine(Real E, Real T, Real Q = 0.0) const override;

  /// helper function to get mu_L from mu_C & neutron energy
  virtual Real getLabCosine(Real E, Real T, Real Q = 0.0) const override;

  /// Function representing the scattering law
  const Function & _scattering_law;
};

#endif // GSL_ENABLED
