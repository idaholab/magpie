/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "DeepNeuralNetFreeEnergy.h"
#include <fstream>

registerADMooseObject("MagpieApp", DeepNeuralNetFreeEnergy);

defineADValidParams(
    DeepNeuralNetFreeEnergy,
    ADMaterial,
    params.addClassDescription(
        "Evaluates a fitted deep neural network to obtain a free energy and its derivatives.");
    params.addParam<FileName>(
        "filename",
        "Data file containing the weights and biasses for a fully connected deep neural network");
    params.addCoupledVar("inputs", "Coupled Variables that are inputs for the neural network");
    params.addParam<std::vector<MaterialPropertyName>>(
        "prop_names", "list of material properties fed from the outputs of the neural network"););

template <ComputeStage compute_stage>
DeepNeuralNetFreeEnergy<compute_stage>::DeepNeuralNetFreeEnergy(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _filename(adGetParam<FileName>("filename")),
    _output_name(adGetParam<std::vector<MaterialPropertyName>>("prop_names")),
    _n_output(_output_name.size()),
    _output(_n_output),
    _n_input(coupledComponents("inputs")),
    _input(_n_input),
    _d_output(_n_input * _n_output)
{
  int x, y;

  // open the NN data file
  std::ifstream ifile;
  ifile.open(_filename);
  if (!ifile)
    paramError("filename", "Unable to open file");

  // read weights (first the number of layers)
  ifile >> _n_layer;
  _weight.resize(_n_layer);
  _z.resize(_n_layer);
  _activation.resize(_n_layer);
  _d_activation.resize(_n_layer);
  for (std::size_t i = 0; i < _n_layer; ++i)
  {
    if (!(ifile >> y >> x))
      mooseError("Error reading file ", _filename);

    // initialize weight matrix
    _weight[i] = DenseMatrix<Real>(x, y);

    // initialize computation buffers (including input and output)
    _activation[i] = DenseVector<ADReal>(y);
    _d_activation[i] = DenseVector<ADReal>(y);
    _z[i] = DenseVector<ADReal>(x);

    for (std::size_t k = 0; k < y; ++k)
      for (std::size_t j = 0; j < x; ++j)
        if (!(ifile >> _weight[i](j, k)))
          mooseError("Error reading weights from file ", _filename);
  }

  // read biases (first the number of layers)
  ifile >> _n_layer;
  _bias.resize(_n_layer);
  for (std::size_t i = 0; i < _n_layer; ++i)
  {
    if (!(ifile >> x))
      mooseError("Error reading file ", _filename);

    _z[i] = DenseVector<ADReal>(x);
    _bias[i] = DenseVector<Real>(x);

    for (std::size_t j = 0; j < x; ++j)
      if (!(ifile >> _bias[i](j)))
        mooseError("Error reading biases from file ", _filename);
  }

  // close parameter file
  ifile.close();

  // validate network properties
  if (_weight.size() != _bias.size())
    paramError("filename", "Inconsistent layer numbers in datafile");
  if (_z[_n_layer - 1].size() != _n_output)
    paramError("prop_names",
               "Number of supplied property names must match the number of output nodes ",
               _z[_n_layer - 1].size(),
               " of the neural net");
  if (_n_input != _activation[0].size())
    paramError("inputs",
               "Number of supplied variables must match the number of input nodes ",
               _activation[0].size(),
               " of the neural net");

  // initialize storage for derivatives
  _prod.resize(_n_layer);
  _diff.resize(_n_layer);
  _diff[0] = _weight[0]; // _prod[0] is not used
  for (std::size_t i = 1; i < _n_layer; ++i)
  {
    _prod[i] = _weight[i];
    _diff[i] = DenseMatrix<ADReal>(_z[i].size(), _n_input);
  }

  // get coupled variables for the inputs
  for (std::size_t j = 0; j < _n_input; ++j)
    _input[j] = &adCoupledValue("inputs", j);

  // create material properties
  std::size_t k = 0;
  for (std::size_t i = 0; i < _n_output; ++i)
  {
    _output[i] = &adDeclareADProperty<Real>(_output_name[i]);
    for (std::size_t j = 0; j < _n_input; ++j)
      _d_output[k++] = &adDeclareADProperty<Real>(
          propertyNameFirst(_output_name[i], this->getVar("inputs", j)->name()));
  }
}

template <ComputeStage compute_stage>
void
DeepNeuralNetFreeEnergy<compute_stage>::computeQpProperties()
{
  // set input nodes
  for (std::size_t j = 0; j < _n_input; ++j)
    _activation[0](j) = (*_input[j])[_qp];

  // evaluate network
  evaluate();

  // copy back output values
  std::size_t k = 0;
  for (std::size_t i = 0; i < _n_output; ++i)
  {
    (*_output[i])[_qp] = _z[_n_layer - 1](i);
    for (std::size_t j = 0; j < _n_input; ++j)
      (*_d_output[k++])[_qp] = _diff[_n_layer - 1](i, j);
  }
}

template <ComputeStage compute_stage>
void
DeepNeuralNetFreeEnergy<compute_stage>::multiply(DenseMatrix<ADReal> & M1,
                                                 const DenseMatrix<ADReal> & M2,
                                                 const DenseMatrix<ADReal> & M3)
{
  // Assertions to make sure we have been
  // passed matrices of the correct dimension.
  libmesh_assert_equal_to(M1.m(), M2.m());
  libmesh_assert_equal_to(M1.n(), M3.n());
  libmesh_assert_equal_to(M2.n(), M3.m());

  const unsigned int m_s = M2.m();
  const unsigned int p_s = M2.n();
  const unsigned int n_s = M1.n();

  for (unsigned int j = 0; j < n_s; j++)
    for (unsigned int i = 0; i < m_s; i++)
    {
      M1.el(i, j) = 0.0;
      for (unsigned int k = 0; k < p_s; k++)
        M1.el(i, j) += M2.el(i, k) * M3.el(k, j);
    }
}

template <ComputeStage compute_stage>
void
DeepNeuralNetFreeEnergy<compute_stage>::evaluate()
{
  std::size_t i = 0;
  while (true)
  {
    // apply weights and biases
    _weight[i].vector_mult(_z[i], _activation[i]);
    _z[i] += _bias[i];

    // derivatives
    if (i > 0)
    {
      // prepare product of weights and activation function derivative (previous layer)
      for (std::size_t j = 0; j < _weight[i].m(); ++j)
        for (std::size_t k = 0; k < _weight[i].n(); ++k)
          _prod[i](j, k) = _weight[i](j, k) * _d_activation[i](k);

      // multiply progressive Jacobian
      multiply(_diff[i], _prod[i], _diff[i - 1]);
    }

    // bail to avoid applying sigmoid to the output
    if (i + 1 == _n_layer)
      break;

    // apply sigmoid activation function
    for (std::size_t j = 0; j < _z[i].size(); ++j)
    {
      _activation[i + 1](j) = 1.0 / (1.0 + std::exp(-_z[i](j)));

      // Note ds(x)/dx = s(x)*(1-s(x))
      // the expensive sigmoid only has to be computed once!
      _d_activation[i + 1](j) = _activation[i + 1](j) * (1 - _activation[i + 1](j));
    }

    // next layer
    ++i;
  }
}
