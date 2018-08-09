/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef INELASTICRECOIL_H
#define INELASTICRECOIL_H

#include "ScatteringRecoilCrossSection.h"

class InelasticRecoil;

template <>
InputParameters validParams<InelasticRecoil>();

class InelasticRecoil : public ScatteringRecoilCrossSection
{
public:
  InelasticRecoil(const InputParameters & parameters);

  virtual void execute() override;

protected:
  /// helper function to get mu_C from T and E
  virtual Real getCMCosine(Real E, Real T, Real Q = 0.0) const override;

  /// helper function to get mu_L from mu_C & neutron energy
  virtual Real getLabCosine(Real E, Real T, Real Q = 0.0) const override;

  /// helper function for maximum recoil energy
  Real getMaxRecoilEnergy(Real E, Real Et);

  /// helper function for maximum recoil energy
  Real getMinRecoilEnergy(Real E, Real Et);

  /// the q values for
  const std::vector<Real> _q_values;

  /// the number of excitation levels == _q_values.size() == _scattering_cross_section.size()
  const unsigned int _n_levels;
};

#endif
