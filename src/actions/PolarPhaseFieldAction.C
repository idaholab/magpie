/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "PolarPhaseFieldAction.h"
#include "FEProblemBase.h"
#include "Factory.h"

registerMooseAction("MagpieApp", PolarPhaseFieldAction, "add_kernel");
registerMooseAction("MagpieApp", PolarPhaseFieldAction, "add_material");

InputParameters
PolarPhaseFieldAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addParam<VariableName>("theta", "theta", "Theta order parameter");
  params.addParam<VariableName>("upsilon", "Upsilon", "Upsilon order parameter");

  params.addParam<Real>("L_theta", 1.0, "Mobility for the theta order parameter");
  params.addParam<Real>("L_upsilon", 1.0, "Mobility for the upsilon order parameter");

  params.addRequiredParam<Real>("a0", "Interpolation coefficient a0");
  params.addRequiredParam<Real>("a_A", "Interpolation coefficient a_A");

  params.addRequiredParam<Real>("a_beta", "Interpolation coefficient a_beta");
  params.addRequiredParam<Real>("a_theta", "Interpolation coefficient a_theta");
  params.addRequiredParam<Real>("a_phi", "Interpolation coefficient a_phi");

  params.addRequiredParam<Real>("beta10", "Gradient energy coefficient between solid 1 and melt");
  params.addRequiredParam<Real>("beta20", "Gradient energy coefficient between solid 2 and melt");
  params.addRequiredParam<Real>("beta21",
                                "Gradient energy coefficient between solid 2 and solid 1");

  params.addRequiredParam<Real>("A10", "Barrier coefficient solid 1 and melt");
  params.addRequiredParam<Real>("A20", "Barrier coefficient solid 2 and melt");
  params.addRequiredParam<Real>("A21", "Barrier coefficient solid 2 and solid 1");

  params.addRequiredParam<Real>("G0", "Thermal energy of the melt");
  params.addRequiredParam<Real>("DeltaG10",
                                "Difference in thermal energy between solid 1 and melt");
  params.addRequiredParam<Real>("DeltaG20",
                                "Difference in thermal energy between solid 2 and melt");

  return params;
}

PolarPhaseFieldAction::PolarPhaseFieldAction(const InputParameters & parameters)
  : Action(parameters)
{
}

void
PolarPhaseFieldAction::act()
{
  // material property names
  const std::string psiL = _name + "_psiL";
  const std::string betaS0 = _name + "_betaS0";
  const std::string beta21phi = _name + "_beta21phi";

  // variables
  auto theta = getParam<VariableName>("theta");
  auto upsilon = getParam<VariableName>("upsilon");

  if (_current_task == "add_material")
  {
    const std::string prefix = _name + "_";
    {
      std::string name = "PolarPFMBetaS0";
      auto params = _factory.getValidParams(name);
      params.set<std::string>("f_name") = betaS0;
      params.set<std::vector<VariableName>>("theta") = {theta};
      params.applyParameters(parameters());
      _problem->addMaterial(name, _name + "_" + name, params);
    }
    {
      std::string name = "PolarPFMPhi";
      auto params = _factory.getValidParams(name);
      params.set<std::string>("f_name") = beta21phi;
      params.set<std::vector<VariableName>>("upsilon") = {upsilon};
      params.applyParameters(parameters());
      _problem->addMaterial(name, _name + "_" + name, params);
    }
    {
      std::string name = "PolarPFMPsiL";
      auto params = _factory.getValidParams(name);
      params.set<std::string>("f_name") = psiL;
      params.set<std::vector<VariableName>>("theta") = {theta};
      params.set<std::vector<VariableName>>("upsilon") = {upsilon};
      params.applyParameters(parameters());
      _problem->addMaterial(name, _name + "_" + name, params);
    }
  }

  if (_current_task == "add_kernel")
  {
    // upsilon
    {
      const std::string prefix = _name + "_upsilon_";
      {
        std::string name = "CoefTimeDerivative";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = upsilon;
        params.set<Real>("Coefficient") = 1.0 / getParam<Real>("L_upsilon");
        _problem->addKernel(name, prefix + name, params);
      }
      {
        std::string name = "PolarPFMDerivative";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = upsilon;
        params.set<MaterialPropertyName>("F") = psiL;
        _problem->addKernel(name, prefix + name, params);
      }
      {
        std::string name = "PolarPFMGradient";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = upsilon;
        params.set<std::vector<VariableName>>("v") = {theta};
        params.set<MaterialPropertyName>("F") = beta21phi;
        _problem->addKernel(name, prefix + name, params);
      }
      {
        std::string name = "MatDiffusion";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = upsilon;
        params.set<MaterialPropertyName>("diffusivity") = betaS0;
        params.set<std::vector<VariableName>>("args") = {theta};
        _problem->addKernel(name, prefix + name, params);
      }
    }

    // theta
    {
      const std::string prefix = _name + "_theta_";
      {
        std::string name = "CoefTimeDerivative";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = theta;
        params.set<Real>("Coefficient") = 1.0 / getParam<Real>("L_theta");
        _problem->addKernel(name, prefix + name, params);
      }
      {
        std::string name = "PolarPFMDerivative";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = theta;
        params.set<MaterialPropertyName>("F") = psiL;
        _problem->addKernel(name, prefix + name, params);
      }
      {
        std::string name = "PolarPFMGradient";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = theta;
        params.set<std::vector<VariableName>>("v") = {theta};
        params.set<MaterialPropertyName>("F") = betaS0;
        _problem->addKernel(name, prefix + name, params);
      }
      {
        std::string name = "MatDiffusion";
        auto params = _factory.getValidParams(name);
        params.set<NonlinearVariableName>("variable") = theta;
        params.set<MaterialPropertyName>("diffusivity") = beta21phi;
        params.set<std::vector<VariableName>>("args") = {upsilon};
        _problem->addKernel(name, prefix + name, params);
      }
    }
  }
}
