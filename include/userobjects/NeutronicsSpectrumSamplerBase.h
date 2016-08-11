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
#ifndef NEUTRONICSSPECTRUMSAMPLERBASE_H
#define NEUTRONICSSPECTRUMSAMPLERBASE_H

#ifdef RATTLESNAKE_ENABLED
  #include "AQData.h"
#endif

// MOOSE includes
#include "ElementUserObject.h"
#include "MultiIndex.h"

// Forward Declarations
class NeutronicsSpectrumSamplerBase;

template<>
InputParameters validParams<NeutronicsSpectrumSamplerBase>();

/**
 * Computes PDFs from neutronics calculations that are
 * used to sample PKAs that will be passed to BCMC simulations.
 */
class NeutronicsSpectrumSamplerBase : public ElementUserObject
{
public:
  NeutronicsSpectrumSamplerBase(const InputParameters & parameters);

  virtual void execute();
  virtual void initialSetup();
  virtual void initialize();
  virtual void finalize();
  virtual void threadJoin(const UserObject & y);
  virtual void meshChanged();

  /// returns a MultiIndex<Real> PDF at a given point ID
  virtual MultiIndex<Real> getPDF(unsigned int point_id) const;

  /// returns a MultiIndex<Real> PDF at a given point ID
  virtual MultiIndex<Real> getPDF(unsigned int point_id) const;
  /// returns a Real magnitude at a given point ID
  virtual Real getMagnitude(unsigned int point_id) const;
  /// returns a std::vector<unsigned int> of ZAIDs at a given point ID
  virtual std::vector<unsigned int> getZAIDs(unsigned int point_id) const;
  /// returns a std::vector<Real> of energies at a given point ID
  virtual std::vector<Real> getEnergies(unsigned int point_id) const;

protected:
  /// a callback executed right before computeRadiationDamagePDF
  virtual void preComputeRadiationDamagePDF();

  /// computes the PKA for isotope i, group g, and SH indices l, m
  virtual Real computeRadiationDamagePDF(unsigned int i, unsigned int g, unsigned int p) = 0;

  /// a subsitute to convert isotope names to zaid if RSN is not available
  unsigned int localStringToZaid(std::string s) const;

  /// vector of target zaids
  const std::vector<std::string> & _target_isotope_names;

  /// the number densities of these isotopes given as variables
  std::vector<const VariableValue *> _number_densities;
  const std::vector<Real> & _energy_group_boundaries;

  /// number of isotopes
  unsigned int _I;
  /// number of energy groups
  unsigned int _G;
  /// spherical harmonics order
  unsigned int _L;

  /// total number of spherical harmonics. Depends on dim.
  unsigned int _nSH;

  /// the points at which PDFs are computed
  const std::vector<Point> & _points;

  /// number of points
  unsigned int _npoints;

  /// flag indicating of recaching the _qp is necessary
  bool _qp_is_cached;

  /// array storing the element id for each point
  std::vector<dof_id_type> _point_element;

  /// the array stores the _qp index for each point
  std::vector<unsigned int> _qp_cache;

  /// stores the radiation damage PDF
  std::vector<MultiIndex<Real> > _sample_point_data;

  /// the current quadrature point
  unsigned int _qp;

  /// the current point
  unsigned int _current_point;

  /// vector of magnitudes at each point ID
  std::vector<Real> _magnitude;
  /// vector of a vector of ZAIDs at each point ID
  std::vector<std::vector<unsigned int> > _zaids;
  /// vector of a vector of Energies at each point ID
  std::vector<std::vector<Real> > _energies;
};

#endif //NEUTRONICSSPECTRUMSAMPLERBASE_H
