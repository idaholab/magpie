/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef RADIALGREENSSOURCE_H
#define RADIALGREENSSOURCE_H

#include "Kernel.h"
#include "RadialGreensConvolution.h"

class RadialGreensSource;

template <>
InputParameters validParams<RadialGreensSource>();

/**
 * Apply the convolution from a RadialGreensConvolution object to a non-linear variable
 */
class RadialGreensSource : public Kernel
{
public:
  RadialGreensSource(const InputParameters & parameters);

protected:
  void precalculateResidual() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  virtual Real computeQpCJacobian();

  // convolution result
  const RadialGreensConvolution::Result & _convolution;

  /// is the kernel used in a coupled form?
  const bool _is_coupled;

  /// int label for the Concentration
  unsigned int _c_var;

  /// Variable value for the concentration
  const VariableValue & _c;

  // rate factor
  const Real _gamma;

  // iterator pointing to the map entry for the current element
  RadialGreensConvolution::Result::const_iterator _result;

  // current timestep size
  const Real & _dt;
};

#endif // RADIALGREENSSOURCE_H
