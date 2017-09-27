/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef RATTLESNAKE_ENABLED
#ifndef NEUTRONICSSPECTRUMSAMPLERSN_H
#define NEUTRONICSSPECTRUMSAMPLERSN_H

#include "NeutronicsSpectrumSamplerBase.h"
#include "AQData.h"

// Forward Declarations
class NeutronicsSpectrumSamplerSN;
class ElasticRecoilCrossSectionUserObject;

template<>
InputParameters validParams<NeutronicsSpectrumSamplerSN>();

/**
 * Computes PDFs from neutronics data that is used to sample PKAs due to
 * other reactions (not fission) for coupled BCMC simulations.
 * The reaction creating the PKAs _must_
 * match target_isotope_names with a unique recoil_isotope_names
 */
class NeutronicsSpectrumSamplerSN : public NeutronicsSpectrumSamplerBase
{
public:
  NeutronicsSpectrumSamplerSN(const InputParameters & parameters);

  virtual Real totalRecoilRate(unsigned int point_id, const std::string & target_isotope) const override;

protected:
  /// computes the PDF for isotope i, group g, and SH indices p
  virtual Real computeRadiationDamagePDF(unsigned int i, unsigned int g, unsigned int p, unsigned int q) override;

  /// vector of target zaids
  const std::vector<std::string> & _recoil_isotope_names;
  /// angular quadrature object
  const AngularQuadrature & _aq;
  /// Number of angular directions in _aq
  unsigned int _ndir;
  /// spherical harmonics coefficient object
  const SHCoefficients _shm;
  /// the angular flux
  std::vector<std::vector<const VariableValue *> > _angular_flux;
  /// UserObjects storing the recoil cross sections
  std::vector<std::vector<const ElasticRecoilCrossSectionUserObject *>> _recoil_cross_sections;
};

#endif //NEUTRONICSSPECTRUMSAMPLERSN_H
#endif //RATTLESNAKE_ENABLED
