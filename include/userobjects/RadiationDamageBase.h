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
#ifdef RATTLESNAKE_ENABLED
#ifndef RADIATIONDAMAGEBASE_H
#define RADIATIONDAMAGEBASE_H

#include "AQData.h"

// MOOSE includes
#include "ElementUserObject.h"
#include "MultiIndex.h"

// Forward Declarations
class RadiationDamageBase;

template<>
InputParameters validParams<RadiationDamageBase>();

/**
 * Computes the PKA species/energy/direction distribution
 * at a given set of point.
 */
class RadiationDamageBase : public ElementUserObject
{
public:
  RadiationDamageBase(const InputParameters & parameters);

  virtual void execute();
  virtual void initialize();
  virtual void finalize();
  virtual void threadJoin(const UserObject & y);
  virtual void meshChanged();
  virtual MultiIndex<Real> getPDF(unsigned int point_id) const;
  virtual Real getMagnitude(unsigned int point_id) const;
  virtual std::vector<unsigned int> getZAIDs(unsigned int point_id) const;
  virtual std::vector<Real> getEnergies(unsigned int point_id) const;

protected:

  /// a callback executed right before computePKA
  virtual void preComputePKA();
  /// computes the PKA for isotope i, group g, and SH indices l, m
  virtual Real computePKA(unsigned int i, unsigned int g, unsigned int p) = 0;

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
  /// the points at which PKAs are computed
  const std::vector<Point> & _points;
  /// number of points
  unsigned int _npoints;
  /// flag indicating of recaching the _qp is necessary
  bool _qp_is_cached;
  /// array storing the element id for each point
  std::vector<dof_id_type> _point_element;
  /// the array stores the _qp index for each point
  std::vector<unsigned int> _qp_cache;
  /// stores the PKA distribution
  std::vector<MultiIndex<Real> > _sample_point_data;

  /// the current quadrature point
  unsigned int _qp;
  /// the current point
  unsigned int _current_point;
};

#endif //RADIATIONDAMAGEBASE_H
#endif //RATTLESNAKE_ENABLED
