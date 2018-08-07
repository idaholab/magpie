/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef NEUTRONICSSPECTRUMSAMPLERFISSION_H
#define NEUTRONICSSPECTRUMSAMPLERFISSION_H

#include "NeutronicsSpectrumSamplerBase.h"

// Forward Declarations
class NeutronicsSpectrumSamplerFission;

template <>
InputParameters validParams<NeutronicsSpectrumSamplerFission>();

/**
 * Computes PDFs from neutronics data that is used to sample PKAs due to fission
 * for coupled BCMC simulations.
 * NOTE: Currently fission is assumed to be isotropic in the LAB
 * frame regardless of the incoming energy.
 */
class NeutronicsSpectrumSamplerFission : public NeutronicsSpectrumSamplerBase
{
public:
  NeutronicsSpectrumSamplerFission(const InputParameters & parameters);

  /// returns a MultiIndex<Real> PDF at a given point ID
  virtual MultiIndex<Real> getPDF(unsigned int point_id) const override;

  virtual Real totalRecoilRate(unsigned int point_id,
                               const std::string & target_isotope) const override;

protected:
  /// computes the PDF for isotope i, group g, and SH indices p
  /// NOTE: for fission p does not mateter
  virtual Real computeRadiationDamagePDF(unsigned int i,
                                         unsigned int g,
                                         unsigned int /*p*/,
                                         unsigned int /*q*/) override;

  /// the scalar flux
  std::vector<const VariableValue *> _scalar_flux;
  /// stores the fission cross sections
  std::vector<std::vector<std::vector<Real>>> _fission_cross_section;
};

#endif // NEUTRONCSSPECTRUMSAMPLERFISSION_H
