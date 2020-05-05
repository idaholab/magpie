/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "FFTProblem.h"
#include "FFTBufferBase.h"
#include "MooseFFTVariable.h"

#include "libmesh/system.h"

registerMooseObject("MagpieApp", FFTProblem);

defineLegacyParams(FFTProblem);

InputParameters
FFTProblem::validParams()
{
  InputParameters params = FEProblem::validParams();
  params.addClassDescription("Enhanced FEProblem that supports FFT buffers as variables.");
  params.set<bool>("skip_nl_system_check") = true;
  return params;
}

FFTProblem::FFTProblem(const InputParameters & parameters)
  : FEProblem(parameters), _fft_dummy_system(*this, "FFTSystem")
{
}

FFTProblem::~FFTProblem()
{
  // delete variable objects
  for (auto & var : _fft_vars)
    for (auto & t : var.second)
      delete t;
}

bool
FFTProblem::hasVariable(const std::string & var_name) const
{
  // first check for FFT buffers
  std::vector<UserObject *> objs;
  theWarehouse().query().condition<AttribThread>(0).condition<AttribName>(var_name).queryInto(objs);
  if (!objs.empty() && dynamic_cast<FFTBufferBase<Real> *>(objs[0]))
  {
    mooseInfo("hasVariable returned true for '", var_name, "' for FFTBuffer");
    return true;
  }

  // fall back to regular variable
  if (FEProblem::hasVariable(var_name))
    return true;

  return false;
}

MooseVariableFEBase &
FFTProblem::getVariable(THREAD_ID tid,
                        const std::string & var_name,
                        Moose::VarKindType expected_var_type,
                        Moose::VarFieldType expected_var_field_type)
{
  // first check for FFT buffers
  std::vector<UserObject *> objs;
  theWarehouse().query().condition<AttribThread>(0).condition<AttribName>(var_name).queryInto(objs);
  if (!objs.empty() && dynamic_cast<FFTBufferBase<Real> *>(objs[0]))
  {
    mooseInfo("getVariable is returning a dummy object for '", var_name, "' for FFTBuffer");

    auto & varlist = _fft_vars[var_name];

    // add dummy name into the dummy system
    unsigned int fft_var_number;
    if (varlist.empty())
      fft_var_number = _fft_dummy_system.system().add_variable(var_name, CONSTANT, MONOMIAL);
    else
      fft_var_number = _fft_dummy_system.system().variable_number(var_name);

    if (varlist.size() <= tid)
      varlist.resize(tid + 1);

    auto params = MooseVariableBase::validParams();
    params.set<MooseEnum>("order") = "CONSTANT";
    params.set<MooseEnum>("family") = "MONOMIAL";
    params.set<unsigned int>("_var_num") = fft_var_number;
    params.set<THREAD_ID>("tid") = tid;
    params.set<THREAD_ID>("_tid") = tid;
    params.set<Moose::VarKindType>("_var_kind") = Moose::VarKindType::VAR_AUXILIARY;
    params.set<SystemBase *>("_system_base") = &_fft_dummy_system;
    params.set<MooseApp *>("_moose_app") = &_app;
    params.set<std::string>("_type") = "MooseFFTVariable";
    params.set<std::string>("_object_name") = var_name;
    params.set<FEProblemBase *>("_fe_problem_base") = this;
    varlist[tid] = new MooseFFTVariable(params);

    return *varlist[tid];
  }

  return FEProblem::getVariable(tid, var_name, expected_var_type, expected_var_field_type);
}
