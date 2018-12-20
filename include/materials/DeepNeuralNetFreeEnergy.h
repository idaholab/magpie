/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "DerivativeParsedMaterialHelper.h"

// Forward Declarations
class DeepNeuralNetFreeEnergy;

template <>
InputParameters validParams<DeepNeuralNetFreeEnergy>();

/**
 *
 */
class DeepNeuralNetFreeEnergy : public DerivativeFunctionMaterialBase
{
public:
  DeepNeuralNetFreeEnergy(const InputParameters & parameters);

protected:
  /// evaluate the network using the inputs in _activation[0]
  void evaluate();

  /// matrix multiplication helper
  void multiply(DenseMatrix<Real> & M1, const DenseMatrix<Real> & M2, const DenseMatrix<Real> & M3);

  /// network data file with weights and biases
  const FileName _filename;

  /// number of inputs
  const std::size_t _n_input;

  /// input variables
  std::vector<const VariableValue *> _input;

  /// network weights
  std::vector<DenseMatrix<Real>> _weight;

  /// network weights
  std::vector<DenseVector<Real>> _bias;

  /// number of layers excluding the input layer
  std::size_t _n_layer;

private:
  /// incoming weighted and biased signal (before applying activation function)
  std::vector<DenseVector<Real>> _z;

  /// activation vectors
  std::vector<DenseVector<Real>> _activation;

  /// derivatives of activation vectors
  std::vector<DenseVector<Real>> _d_activation;

  /// progressive jacobian matrices
  std::vector<DenseMatrix<Real>> _diff;

  /// product of the weight matrix and the derivative of the activation function
  std::vector<DenseMatrix<Real>> _prod;
};
