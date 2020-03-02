/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "AuxKernel.h"

class MDRunBase;

class MDGranularPropertyAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MDGranularPropertyAux(const InputParameters & params);

  virtual Real computeValue();

  static MooseEnum mdAveragingType();

protected:
  /// referene to the MDRunBase user object
  const MDRunBase & _md_uo;

  /// the type of average to be computed
  MooseEnum _average_type;

  /// ID of the desired MD property
  unsigned int _property_id;

  /// property value that is computed only on qp = 0
  Real _property_value;
};
