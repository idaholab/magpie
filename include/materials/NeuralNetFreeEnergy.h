/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "NeuralNetFreeEnergyBase.h"

/**
 * Evaluate a deep neural net and its derivatives
 */
class NeuralNetFreeEnergy : public NeuralNetFreeEnergyBase
{
public:
  static InputParameters validParams();

  NeuralNetFreeEnergy(const InputParameters & parameters);

protected:
  enum class ActivationFunction
  {
    SIGMOID,
    SOFTSIGN,
    TANH
  } _activation_function;

  /// apply activation functions (and record their derivatives) for the current layer
  void applyLayerActivation() override;
};
