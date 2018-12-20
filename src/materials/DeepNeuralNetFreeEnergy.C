/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "DeepNeuralNetFreeEnergy.h"
#include <fstream>

registerMooseObject("MagpieApp", DeepNeuralNetFreeEnergy);

template <>
InputParameters
validParams<DeepNeuralNetFreeEnergy>()
{
  InputParameters params = validParams<DerivativeFunctionMaterialBase>();
  params.addClassDescription(
      "Evaluates a fitted deep neural network to obtain a free energy and its derivatives.");
  params.addParam<FileName>(
      "filename",
      "Data file containing the weights and biasses for a fully connected deep neural network");
  params.addCoupledVar("inputs", "Coupled Variables that are inputs for the neural network");
  return params;
}

DeepNeuralNetFreeEnergy::DeepNeuralNetFreeEnergy(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
    _filename(getParam<FileName>("filename")),
    _n_input(coupledComponents("inputs")),
    _input(_n_input)
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
    _activation[i] = DenseVector<Real>(y);
    _d_activation[i] = DenseVector<Real>(y);
    _z[i] = DenseVector<Real>(x);

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

    _z[i] = DenseVector<Real>(x);
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
  if (_z[_n_layer - 1].size() != 1)
    paramError("filename", "Neural network needs to have exactly one output node");
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
    _diff[i] = DenseMatrix<Real>(_z[i].size(), _n_input);
  }

  // test output
  for (Real T = 0.0; T <= 1.0; T += 0.05)
    for (Real c = 0.0; c <= 1.0; c += 0.05)
    {
      _activation[0](0) = T + 0.001;
      _activation[0](1) = c;
      evaluate();
      Real out10 = _z[_n_layer - 1](0);

      _activation[0](0) = T;
      _activation[0](1) = c + 0.001;
      evaluate();
      Real out01 = _z[_n_layer - 1](0);

      _activation[0](0) = T;
      _activation[0](1) = c;
      evaluate();
      Real out00 = _z[_n_layer - 1](0);

      Real fd_T = (out10 - out00) / 0.001;
      Real fd_c = (out01 - out00) / 0.001;

      Real hc_T = _diff[_n_layer - 1](0, 0);
      Real hc_c = _diff[_n_layer - 1](0, 1);

      std::cout << "DNN" << T << ' ' << c << ' ' << out00 << ' ' << hc_T << ' ' << fd_T << ' '
                << hc_c << ' ' << fd_c << '\n';
    }
}

void
DeepNeuralNetFreeEnergy::multiply(DenseMatrix<Real> & M1,
                                  const DenseMatrix<Real> & M2,
                                  const DenseMatrix<Real> & M3)
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

void
DeepNeuralNetFreeEnergy::evaluate()
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
      const std::size_t m = _weight[i].m();
      const std::size_t n = _weight[i].n();
      for (std::size_t j = 0; j < m; ++j)
        for (std::size_t k = 0; k < n; ++k)
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
