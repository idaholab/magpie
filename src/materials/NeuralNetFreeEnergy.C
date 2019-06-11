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

    MooseEnum activationFunctionEnum("SIGMOID", "SIGMOID");
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
        _activation[_layer + 1](j) = 1.0 / (1.0 + std::exp(-_z[_layer](j)));

        // Note ds(x)/dx = s(x)*(1-s(x))
        // the expensive sigmoid only has to be computed once!
        _d_activation[_layer + 1](j) =
            _activation[_layer + 1](j) * (1 - _activation[_layer + 1](j));
      }
      return;

    default:
      paramError("activation_function", "Unknown activation function");
  }
}
