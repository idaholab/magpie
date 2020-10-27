/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "NeuralNetFreeEnergyBase.h"
#include "MooseEnum.h"
#include <fstream>

InputParameters
NeuralNetFreeEnergyBase::validParams()
{
  auto params = ADMaterial::validParams();
  params.addClassDescription(
      "Evaluates a fitted deep neural network to obtain a free energy and its derivatives.");

  MooseEnum fileFormatEnum("MAGPIE GENANN");
  params.template addParam<MooseEnum>(
      "file_format", fileFormatEnum, "Weights and biases file format");

  params.template addParam<FileName>(
      "file_name",
      "Data file containing the weights and biasses for a fully connected deep neural network");
  params.addCoupledVar("inputs", "Coupled Variables that are inputs for the neural network");
  params.template addParam<std::vector<MaterialPropertyName>>(
      "prop_names", "list of material properties fed from the outputs of the neural network");
  params.template addParam<bool>(
      "debug", false, "Tabulate the NN to a file for debugging purposes");
  return params;
}

NeuralNetFreeEnergyBase::NeuralNetFreeEnergyBase(const InputParameters & parameters)
  : ADMaterial(parameters),
    _file_format(getParam<MooseEnum>("file_format").template getEnum<FileFormat>()),
    _file_name(getParam<FileName>("file_name")),
    _output_name(getParam<std::vector<MaterialPropertyName>>("prop_names")),
    _n_output(_output_name.size()),
    _output(_n_output),
    _n_input(coupledComponents("inputs")),
    _input(_n_input),
    _d_output(_n_input * _n_output)
{
  // open the NN data file
  std::ifstream ifile;
  ifile.open(_file_name);
  if (!ifile)
    paramError("file_name", "Unable to open file");

  // call the reader for the requested format
  switch (_file_format)
  {
    case FileFormat::MAGPIE:
      loadMagpieNet(ifile);
      break;

    case FileFormat::GENANN:
      loadGenANN(ifile);
      break;

    default:
      paramError("file_format", "Unknown file format");
  }

  // close parameter file
  ifile.close();

  // validate network properties
  if (_weight.size() != _bias.size())
    paramError("file_name", "Inconsistent layer numbers in datafile");
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
    _output[i] = &declareADProperty<Real>(_output_name[i]);
    for (std::size_t j = 0; j < _n_input; ++j)
      _d_output[k++] = &declareADProperty<Real>(
          derivativePropertyNameFirst(_output_name[i], this->getVar("inputs", j)->name()));
  }
}

void
NeuralNetFreeEnergyBase::initialSetup()
{
  if (getParam<bool>("debug"))
    debugDump();
}

void
NeuralNetFreeEnergyBase::debugDump()
{
  std::ofstream ofile;
  ofile.open(_name + "_debug.dat");
  for (Real T = 0.0; T <= 1.0; T += 0.01)
    for (Real c = 0.0; c <= 1.0; c += 0.01)
    {
      ofile << T << ' ' << c;

      // set inputs
      _activation[0](0) = T;
      _activation[0](1) = c;

      // evaluate network
      evaluate();

      // output values
      for (std::size_t i = 0; i < _n_output; ++i)
        ofile << ' ' << _z[_n_layer - 1](i);

      ofile << '\n';
    }
}

void
NeuralNetFreeEnergyBase::loadGenANN(std::ifstream & ifile)
{
  std::size_t x, y;

  // read layer layout
  unsigned int n_inputs, n_hidden_layers, n_hidden, n_outputs;
  if (!(ifile >> n_inputs >> n_hidden_layers >> n_hidden >> n_outputs))
    mooseError("Error reading genann file header from file ", _file_name);
  _n_layer = n_hidden_layers + 1;
  _bias.resize(_n_layer);
  _weight.resize(_n_layer);
  _z.resize(_n_layer);
  _activation.resize(_n_layer);
  _d_activation.resize(_n_layer);

  // read weights and biases
  for (std::size_t i = 0; i < _n_layer; ++i)
  {
    // negative biases come first
    x = i < _n_layer - 1 ? n_hidden : n_outputs;

    _z[i] = DenseVector<ADReal>(x);
    _bias[i] = DenseVector<Real>(x);
    for (std::size_t j = 0; j < x; ++j)
      if (!(ifile >> _bias[i](j)))
        mooseError("Error reading biases from file ", _file_name);
      else
        _bias[i](j) *= -1.0;

    // next come the weights
    y = i > 0 ? n_hidden : n_inputs;
    // initialize weight matrix
    _weight[i] = DenseMatrix<Real>(x, y);

    // initialize computation buffers (including input and output)
    _activation[i] = DenseVector<ADReal>(y);
    _d_activation[i] = DenseVector<ADReal>(y);
    _z[i] = DenseVector<ADReal>(x);

    for (std::size_t k = 0; k < y; ++k)
      for (std::size_t j = 0; j < x; ++j)
        if (!(ifile >> _weight[i](j, k)))
          mooseError("Error reading weights from file ", _file_name);
  }
}

void
NeuralNetFreeEnergyBase::loadMagpieNet(std::ifstream & ifile)
{
  std::size_t x, y;

  // read weights (first the number of layers)
  ifile >> _n_layer;
  _weight.resize(_n_layer);
  _z.resize(_n_layer);
  _activation.resize(_n_layer);
  _d_activation.resize(_n_layer);
  for (std::size_t i = 0; i < _n_layer; ++i)
  {
    if (!(ifile >> y >> x))
      mooseError("Error reading file ", _file_name);

    // initialize weight matrix
    _weight[i] = DenseMatrix<Real>(x, y);

    // initialize computation buffers (including input and output)
    _activation[i] = DenseVector<ADReal>(y);
    _d_activation[i] = DenseVector<ADReal>(y);
    _z[i] = DenseVector<ADReal>(x);

    for (std::size_t k = 0; k < y; ++k)
      for (std::size_t j = 0; j < x; ++j)
        if (!(ifile >> _weight[i](j, k)))
          mooseError("Error reading weights from file ", _file_name);
  }

  // read biases
  _bias.resize(_n_layer);
  for (std::size_t i = 0; i < _n_layer; ++i)
  {
    if (!(ifile >> x))
      mooseError("Error reading file ", _file_name);

    _z[i] = DenseVector<ADReal>(x);
    _bias[i] = DenseVector<Real>(x);

    for (std::size_t j = 0; j < x; ++j)
      if (!(ifile >> _bias[i](j)))
        mooseError("Error reading biases from file ", _file_name);
  }
}

void
NeuralNetFreeEnergyBase::computeQpProperties()
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

void
NeuralNetFreeEnergyBase::multiply(DenseMatrix<ADReal> & M1,
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

void
NeuralNetFreeEnergyBase::evaluate()
{
  _layer = 0;
  while (true)
  {
    // apply weights and biases
    _weight[_layer].vector_mult(_z[_layer], _activation[_layer]);
    _z[_layer] += _bias[_layer];

    // derivatives
    if (_layer > 0)
    {
      // prepare product of weights and activation function derivative (previous layer)
      for (std::size_t j = 0; j < _weight[_layer].m(); ++j)
        for (std::size_t k = 0; k < _weight[_layer].n(); ++k)
          _prod[_layer](j, k) = _weight[_layer](j, k) * _d_activation[_layer](k);

      // multiply progressive Jacobian
      multiply(_diff[_layer], _prod[_layer], _diff[_layer - 1]);
    }

    // bail to avoid applying activation function to the output
    if (_layer + 1 == _n_layer)
      break;

    // apply activation function
    applyLayerActivation();

    // next layer
    ++_layer;
  }
}
