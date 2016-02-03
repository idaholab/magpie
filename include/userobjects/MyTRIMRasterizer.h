/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MYTRIMRASTERIZER_H
#define MYTRIMRASTERIZER_H

#include "ElementUserObject.h"
#include "PKAGeneratorBase.h"

#include <map>
#include <vector>

class MyTRIMRasterizer;

template<>
InputParameters validParams<MyTRIMRasterizer>();

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library
 */
class MyTRIMRasterizer : public ElementUserObject
{
public:
  MyTRIMRasterizer(const InputParameters & parameters);

  /// determines if a TRIM run is executed during this timestep
  virtual bool executeThisTimestep() const;

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual void finalize();

  // get the concentration array
  const std::vector<Real> & material(const Elem *) const;

  // get the mass array
  const std::vector<Real> & mass() const { return _trim_mass; }

  // get the charge array
  const std::vector<Real> & charge() const { return _trim_charge; }

  // get the PKA list
  const std::vector<MyTRIM_NS::IonBase> & getPKAList() const { return _pka_list; }

  // get the variable ID of the first coupled variable (to determine the periodicity)
  int periodic() const { return _periodic; }

  // get the number of elements in the TRIM simulation
  unsigned int nVars() const { return _nvars; }

protected:
  /// number of coupled variables to map
  const unsigned int _nvars;

  ///@{ Element data
  std::vector<Real> _trim_mass;
  std::vector<Real> _trim_charge;
  ///@}

  /// coupled variable values
  std::vector<const VariableValue *> _var;

  /// @{ PKA generators
  const std::vector<UserObjectName> _pka_generator_names;
  std::vector<const PKAGeneratorBase *> _pka_generators;
  /// @}

  /// material map for the TRIM simulation
  typedef std::map<dof_id_type, std::vector<Real> > MaterialMap;
  MaterialMap _material_map;

  /// variable number to use for minPeriodicDistance calls (i.e. use the periodicity of this variable)
  const int _periodic;

  /// cumulative PKA list
  std::vector<MyTRIM_NS::IonBase> _pka_list;

  /// last time the BCMC simulation ran
  Real _last_time;

  /// End time of teh curent step
  Real _step_end_time;

private:
  bool _execute_this_timestep;
};

#endif //MYTRIMRASTERIZER_H
