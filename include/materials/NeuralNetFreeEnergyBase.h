/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "ADMaterial.h"
#include "DerivativeMaterialPropertyNameInterface.h"

#define usingNeuralNetFreeEnergyBaseMembers                                                        \
  usingMaterialMembers;                                                                            \
  using NeuralNetFreeEnergyBase<compute_stage>::_layer;                                            \
  using NeuralNetFreeEnergyBase<compute_stage>::_z;                                                \
  using NeuralNetFreeEnergyBase<compute_stage>::_activation;                                       \
  using NeuralNetFreeEnergyBase<compute_stage>::_d_activation;                                     \
  using NeuralNetFreeEnergyBase<compute_stage>::applyLayerActivation

// Forward Declarations
template <ComputeStage>
class NeuralNetFreeEnergyBase;

declareADValidParams(NeuralNetFreeEnergyBase);

/**
 * Evaluate a deep neural net and its derivatives
 */
template <ComputeStage compute_stage>
class NeuralNetFreeEnergyBase : public ADMaterial<compute_stage>,
                                public DerivativeMaterialPropertyNameInterface
{
public:
  NeuralNetFreeEnergyBase(const InputParameters & parameters);

  virtual void initialSetup();

protected:
  /// compute material properties for the current quadrature point
  virtual void computeQpProperties();

  /// apply activation functions (and record their derivatives) for the current layer
  virtual void applyLayerActivation() = 0;

  /// layer currently processed
  std::size_t _layer;

  /// incoming weighted and biased signal (before applying activation function)
  std::vector<DenseVector<ADReal>> _z;

  /// activation vectors
  std::vector<DenseVector<ADReal>> _activation;

  /// derivatives of activation vectors
  std::vector<DenseVector<ADReal>> _d_activation;

private:
  /// load a neural net description form a genann file
  void loadGenANN(std::ifstream & ifile);

  /// load a neural net description form a custom weights and biases file
  void loadMagpieNet(std::ifstream & ifile);

  /// matrix multiplication helper
  void multiply(DenseMatrix<ADReal> & M1,
                const DenseMatrix<ADReal> & M2,
                const DenseMatrix<ADReal> & M3);

  /// debug output of the NN results
  void debugDump();

  /// evaluate the network using the inputs in _activation[0]
  void evaluate();

  /// format of the file containing the weights and biases data
  enum class FileFormat
  {
    MAGPIE,
    GENANN
  } _file_format;

  /// network data file with weights and biases
  const FileName _file_name;

  /// output material property names
  const std::vector<MaterialPropertyName> _output_name;

  /// number of outputs
  const std::size_t _n_output;

  /// output properties
  std::vector<ADMaterialProperty(Real) *> _output;

  /// number of inputs
  const std::size_t _n_input;

  /// network weights
  std::vector<DenseMatrix<Real>> _weight;

  /// network weights
  std::vector<DenseVector<Real>> _bias;

  /// number of layers excluding the input layer
  std::size_t _n_layer;

  /// input variables
  std::vector<const ADVariableValue *> _input;

  /// output property derivatives
  std::vector<ADMaterialProperty(Real) *> _d_output;

  /// progressive Jacobian matrices
  std::vector<DenseMatrix<ADReal>> _diff;

  /// product of the weight matrix and the derivative of the activation function
  std::vector<DenseMatrix<ADReal>> _prod;

  usingMaterialMembers;
};
