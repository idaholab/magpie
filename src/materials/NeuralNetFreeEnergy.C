/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "NeuralNetFreeEnergy.h"

registerADMooseObject("MagpieApp", NeuralNetFreeEnergy);

defineADValidParams(
    NeuralNetFreeEnergy,
    NeuralNetFreeEnergyBase,
    params.addClassDescription("Evaluates a fitted deep neural network to obtain a free energy and "
                               "its derivatives with a preset activation function.");

    MooseEnum activationFunctionEnum("SIGMOID SOFTSIGN TANH", "SIGMOID");
    params.addParam<MooseEnum>("activation_function",
                               activationFunctionEnum,
                               "Weights and biases file format"););

template <ComputeStage compute_stage>
NeuralNetFreeEnergy<compute_stage>::NeuralNetFreeEnergy(const InputParameters & parameters)
  : NeuralNetFreeEnergyBase<compute_stage>(parameters),
    _activation_function(
        getParam<MooseEnum>("activation_function").template getEnum<ActivationFunction>())
{
}

template <ComputeStage compute_stage>
void
NeuralNetFreeEnergy<compute_stage>::applyLayerActivation()
{
  switch (_activation_function)
  {
    case ActivationFunction::SIGMOID:
      for (std::size_t j = 0; j < _z[_layer].size(); ++j)
      {
        const auto & z = _z[_layer](j);

        const auto F = 1.0 / (1.0 + std::exp(-z));
        _activation[_layer + 1](j) = F;

        // Note dF(z)/dz = F(z)*(1-F(z)), thus the expensive sigmoid only has to be computed once!
        _d_activation[_layer + 1](j) = F * (1 - F);
      }
      return;

    case ActivationFunction::SOFTSIGN:
      for (std::size_t j = 0; j < _z[_layer].size(); ++j)
      {
        const auto & z = _z[_layer](j);

        const auto p = 1.0 + std::abs(z);
        const auto F = z / p;
        _activation[_layer + 1](j) = F;

        const auto dF = -std::abs(z) / (p * p) + 1.0 / p;
        _d_activation[_layer + 1](j) = dF;
      }
      return;

    case ActivationFunction::TANH:
      for (std::size_t j = 0; j < _z[_layer].size(); ++j)
      {
        const auto & z = _z[_layer](j);

        const auto F = std::tanh(z);
        _activation[_layer + 1](j) = F;

        _d_activation[_layer + 1](j) = 1.0 - F * F;
      }
      return;

    default:
      paramError("activation_function", "Unknown activation function");
  }
}
