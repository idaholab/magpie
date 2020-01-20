/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "FEProblem.h"
#include "AuxiliarySystem.h"

/**
 * Enhanced FEProblem that supports FFT buffers as variables
 */
class FFTProblem : public FEProblem
{
public:
  static InputParameters validParams();

  FFTProblem(const InputParameters & parameters);
  ~FFTProblem();

  virtual bool hasVariable(const std::string & var_name) const override;
  virtual MooseVariableFEBase & getVariable(
      THREAD_ID tid,
      const std::string & var_name,
      Moose::VarKindType expected_var_type = Moose::VarKindType::VAR_ANY,
      Moose::VarFieldType expected_var_field_type = Moose::VarFieldType::VAR_FIELD_ANY) override;

  // getMaterialData

protected:
  /// map from variable name to list of variable objects (one per thread)
  std::map<std::string, std::vector<MooseVariableFEBase *>> _fft_vars;

  /// dummy system for the FFT variables
  AuxiliarySystem _fft_dummy_system;

  unsigned int _fft_var_number;
};
